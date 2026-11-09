import os
import requests
import time
import random
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Set
from urllib.parse import quote_plus
import re
from datetime import datetime
import json
from collections import defaultdict
import pandas as pd
from elsapy.elsclient import ElsClient
from elsapy.elssearch import ElsSearch
import traceback
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm


class RateLimitedSession:
    """Session with rate limiting and retry functionality."""
    def __init__(self, calls_per_second: float = 1):
        self.session = requests.Session()
        self.calls_per_second = calls_per_second
        self.last_call = 0
        retry_strategy = Retry(
            total=3, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)

    def request(self, *args, **kwargs):
        current_time = time.time()
        time_since_last_call = current_time - self.last_call
        wait_time = (1 / self.calls_per_second) - time_since_last_call
        if wait_time > 0:
            time.sleep(wait_time)

        response = self.session.request(*args, **kwargs)
        self.last_call = time.time()
        return response


class PubChemDB:
    """Interface for PubChem database."""
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.session = RateLimitedSession(calls_per_second=4.0)

    def download_3d_structure(self, cid: str) -> Optional[bytes]:
        try:
            url = f"{self.base_url}/compound/cid/{cid}/SDF?record_type=3d"
            response = self.session.request("GET", url, timeout=30)
            response.raise_for_status()
            if response.content and b"$$$$" in response.content:
                return response.content
            else:
                logging.warning(f"PubChem download for CID {cid} did not return valid SDF content.")
                return None
        except requests.exceptions.RequestException as e:
            status_code = e.response.status_code if e.response is not None else 'N/A'
            logging.warning(f"PubChem download HTTP error for CID {cid} (Status: {status_code}): {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected PubChem download error for CID {cid}: {e}")
            return None


class SDFDownloader:
    """Downloads 3D SDF structures from PubChem using CIDs."""
    def __init__(self, max_workers: int = 4):
        self.max_workers = max_workers
        self._setup_logging()
        self.database = PubChemDB()
        self.downloaded_files: Dict[str, Path] = {}
        self.failed_cids: Set[str] = set()

    def _setup_logging(self):
        log_dir = Path("logs")
        log_dir.mkdir(exist_ok=True)
        self.log_file = log_dir / 'sdf_download.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s [%(levelname)s] (%(threadName)s) %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler()
            ]
        )
        logging.info("--- Starting new SDF download session ---")

    def save_structure(self, content: bytes, cid: str, output_dir: Path) -> Optional[Path]:
        """Saves structure content to a file."""
        try:
            filename = f"{cid}.sdf"
            filepath = output_dir / filename
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(filepath, 'wb') as f:
                f.write(content)
            logging.info(f"Saved SDF for CID {cid} to {filepath}")
            return filepath
        except IOError as e:
            logging.error(f"I/O error saving SDF for CID {cid}: {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected error saving SDF for CID {cid}: {e}")
            return None

    def download_sdf(self, cid: str, output_dir: Path) -> Optional[Path]:
        """Downloads SDF for a single CID."""
        logging.debug(f"Downloading SDF for CID: {cid}")
        content = self.database.download_3d_structure(cid)
        if content:
            saved_path = self.save_structure(content, cid, output_dir)
            if saved_path:
                return saved_path
        self.failed_cids.add(cid)
        return None

    def process_cid_list(self, cids: List[str], output_dir_name: str):
        """Processes a list of CIDs, downloading SDFs to a specific directory."""
        specific_output_dir = Path("02_Ligands") / "raw_sdf" / output_dir_name
        specific_output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"--- Processing CIDs for output directory: {specific_output_dir} ---")

        self.downloaded_files.clear()
        self.failed_cids.clear()

        with ThreadPoolExecutor(max_workers=self.max_workers, thread_name_prefix='Downloader') as executor:
            future_to_cid = {
                executor.submit(self.download_sdf, cid, specific_output_dir): cid
                for cid in set(cids)
            }

            total_cids = len(future_to_cid)
            with tqdm(total=total_cids, desc=f"Downloading for {output_dir_name}", unit="cid") as pbar:
                for future in as_completed(future_to_cid):
                    cid = future_to_cid[future]
                    try:
                        result_path = future.result()
                        if result_path:
                            self.downloaded_files[cid] = result_path
                    except Exception as e:
                        logging.error(f"Thread error processing CID {cid}: {e}")
                        self.failed_cids.add(cid)
                    finally:
                        pbar.update(1)

        logging.info(f"--- Finished processing for {specific_output_dir} ---")
        successfully_downloaded = len(self.downloaded_files)
        failed_count = len(self.failed_cids)

        logging.info(f"Total unique CIDs attempted: {total_cids}")
        logging.info(f"Successfully downloaded: {successfully_downloaded}")
        logging.info(f"Failed to download: {failed_count}")

        failed_list_filename = specific_output_dir / "failed_downloads.txt"
        try:
            with open(failed_list_filename, 'w', encoding='utf-8') as f:
                if self.failed_cids:
                    f.write("# CIDs that could not be downloaded\n")
                    for cid in sorted(list(self.failed_cids)):
                        f.write(f"{cid}\n")
                else:
                    f.write("# No CIDs failed to download.\n")
            logging.info(f"List of failed CIDs saved to: {failed_list_filename}")
        except IOError as e:
            logging.error(f"Could not write failed CIDs list {failed_list_filename}: {e}")

        return {
            "total_attempted": total_cids,
            "downloaded": successfully_downloaded,
            "failed": failed_count,
            "output_dir": str(specific_output_dir),
            "failed_list_file": str(failed_list_filename)
        }


def main():
    """Main function to process the Excel files and download SDFs."""
    excel_files = [
        {
            "file": "01_Literature_Data/manual_curation/endophyte_derived_for_d_anomala.csv.xlsx",
            "output_dir": "endophyte_derived"
        },
        {
            "file": "01_Literature_Data/manual_curation/high_confidence_compounds.csv.xlsx",
            "output_dir": "high_confidence"
        },
        {
            "file": "01_Literature_Data/manual_curation/require_structural_curation.csv.xlsx",
            "output_dir": "require_structural_curation"
        }
    ]

    downloader = SDFDownloader(max_workers=4)
    all_run_stats = []

    for config in excel_files:
        file_path = config["file"]
        output_dir = config["output_dir"]

        print(f"\n===== Processing: {file_path} =====")
        if not Path(file_path).exists():
            print(f"Error: Excel file not found: {file_path}")
            logging.error(f"Excel file not found: {file_path}, skipping.")
            all_run_stats.append({**config, "status": "skipped - file not found"})
            continue

        try:
            df = pd.read_excel(file_path)
            if 'cid' not in df.columns:
                print(f"Error: 'cid' column not found in {file_path}")
                logging.error(f"'cid' column not found in {file_path}, skipping.")
                all_run_stats.append({**config, "status": "skipped - no cid column"})
                continue

            cids = df['cid'].dropna().astype(str).str.strip().tolist()
            unique_cids = list(set(cids))
            if not unique_cids:
                print(f"Warning: No valid CIDs found in {file_path}")
                logging.warning(f"No valid CIDs in {file_path}, skipping.")
                all_run_stats.append({**config, "status": "skipped - no cids"})
                continue

            print(f"Found {len(unique_cids)} unique CIDs in {file_path}.")

            run_stats = downloader.process_cid_list(unique_cids, output_dir)
            all_run_stats.append({**config, **run_stats, "status": "completed"})

        except Exception as e:
            print(f"An error occurred while processing {file_path}: {e}")
            logging.exception(f"Error processing {file_path}")
            all_run_stats.append({**config, "status": f"error - {e}"})

    print("\n\n===== Overall SDF Download Summary =====")
    total_attempted_all = 0
    total_downloaded_all = 0
    total_failed_all = 0

    for stats in all_run_stats:
        print(f"\n--- Summary for: {stats['file']} ---")
        print(f"  Status: {stats['status']}")
        if stats.get('status') == 'completed':
            print(f"  Output Directory: {stats['output_dir']}")
            print(f"  CIDs Attempted: {stats['total_attempted']}")
            print(f"  Successfully Downloaded: {stats['downloaded']}")
            print(f"  Failed to Download: {stats['failed']}")
            print(f"  Failed List File: {stats['failed_list_file']}")
            total_attempted_all += stats['total_attempted']
            total_downloaded_all += stats['downloaded']
            total_failed_all += stats['failed']
        else:
            print(f"  Processing was skipped or errored.")

    print("\n--- Grand Totals ---")
    print(f"Total unique CIDs attempted across all files: {total_attempted_all}")
    print(f"Total SDFs successfully downloaded across all files: {total_downloaded_all}")
    print(f"Total unique CIDs failed across all files: {total_failed_all}")
    print(f"\nCheck '{downloader.log_file}' for detailed logging information.")
    print("===================================")


if __name__ == "__main__":
    main()
