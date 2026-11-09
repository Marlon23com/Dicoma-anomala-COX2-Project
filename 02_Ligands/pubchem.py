import requests
import os
import time
import logging
import re
from pathlib import Path
from typing import Optional, Dict, List, Set
import json
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from urllib.parse import quote
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

class ChemicalDatabase(ABC):
    """Abstract base class for chemical database interfaces."""
    @abstractmethod
    def search_compound(self, name: str) -> Optional[str]:
        pass

    @abstractmethod
    def download_3d_structure(self, db_id: str, name: str) -> Optional[bytes]:
        pass

class PubChemDB(ChemicalDatabase):
    """Interface for PubChem database."""
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.session = RateLimitedSession(calls_per_second=4.0)

    def search_compound(self, name: str) -> Optional[str]:
        try:
            url = f"{self.base_url}/compound/name/{quote(name)}/cids/JSON"
            response = self.session.request("GET", url, timeout=20)
            response.raise_for_status()
            data = response.json()
            if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
                return str(data['IdentifierList']['CID'][0])
            logging.debug(f"No CID found for {name} in PubChem.")
            return None
        except requests.exceptions.RequestException as e:
            # Specifically log HTTP errors like 404 Not Found
            status_code = e.response.status_code if e.response is not None else 'N/A'
            logging.warning(f"PubChem search HTTP error for '{name}' (Status: {status_code}): {e}")
            return None
        except json.JSONDecodeError as e:
            logging.error(f"PubChem search JSON decode error for '{name}': {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected PubChem search error for '{name}': {e}")
            return None

    def download_3d_structure(self, db_id: str, name: str) -> Optional[bytes]:
        try:
            url = f"{self.base_url}/compound/cid/{db_id}/SDF?record_type=3d"
            response = self.session.request("GET", url, timeout=30)
            response.raise_for_status()
            if response.content and b"$$$$" in response.content:
                 return response.content
            else:
                 logging.warning(f"PubChem download for {name} (CID: {db_id}) did not return valid SDF content.")
                 return None
        except requests.exceptions.RequestException as e:
            status_code = e.response.status_code if e.response is not None else 'N/A'
            logging.warning(f"PubChem download HTTP error for {name} (CID: {db_id}, Status: {status_code}): {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected PubChem download error for {name} (CID: {db_id}): {e}")
            return None

class ChemSpiderDB(ChemicalDatabase):
    """Interface for ChemSpider database."""
    def __init__(self, api_key: str):
        self.base_url = "https://api.rsc.org/compounds/v1"
        # Check ChemSpider rate limits, adjust if needed
        self.session = RateLimitedSession(calls_per_second=0.5)
        if not api_key:
            logging.warning("ChemSpider API key not provided. ChemSpider functionality will be limited.")
            self.api_key = None
            self.headers = {'Accept': 'application/json'} # Basic header without API key
        else:
            self.api_key = api_key
            self.headers = {'apikey': self.api_key, 'Accept': 'application/json'}

    def search_compound(self, name: str) -> Optional[str]:
        if not self.api_key: return None # Skip if no API key
        try:
            # ChemSpider uses a POST request for batch search, even for one item
            search_url = f"{self.base_url}/filter/name"
            payload = {"name": name, "orderBy": "", "orderDirection": ""} # Simplified payload
            response = self.session.request(
                "POST", search_url, headers=self.headers, json=payload, timeout=20
            )
            response.raise_for_status()
            data = response.json()
            # The queryId is needed for the next step to get status/results
            query_id = data.get('queryId')
            if not query_id:
                 logging.debug(f"No queryId returned from ChemSpider search for {name}")
                 return None

            status_url = f"{self.base_url}/filter/{query_id}/status"
            for _ in range(5): # Poll status a few times
                time.sleep(2) # Wait for search to potentially complete
                status_response = self.session.request("GET", status_url, headers=self.headers, timeout=15)
                status_response.raise_for_status()
                status_data = status_response.json()
                status = status_data.get('status')
                if status == 'Complete':
                    # Get results
                    results_url = f"{self.base_url}/filter/{query_id}/results"
                    results_response = self.session.request("GET", results_url, headers=self.headers, timeout=15)
                    results_response.raise_for_status()
                    results_data = results_response.json()
                    if results_data.get('results'):
                        # Return the first result ID (assuming it's the best match)
                        return str(results_data['results'][0])
                    else:
                        logging.debug(f"No results found in ChemSpider for queryId {query_id} ({name})")
                        return None
                elif status == 'Failed' or status == 'Suspended':
                    logging.warning(f"ChemSpider search failed or suspended for {name} (QueryId: {query_id}, Status: {status})")
                    return None
            logging.warning(f"ChemSpider search timed out waiting for status update for {name} (QueryId: {query_id})")
            return None

        except requests.exceptions.RequestException as e:
            status_code = e.response.status_code if e.response is not None else 'N/A'
            logging.warning(f"ChemSpider search HTTP error for '{name}' (Status: {status_code}): {e}")
            return None
        except json.JSONDecodeError as e:
            logging.error(f"ChemSpider search JSON decode error for '{name}': {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected ChemSpider search error for '{name}': {e}")
            return None

    def download_3d_structure(self, db_id: str, name: str) -> Optional[bytes]:
        if not self.api_key: return None # Skip if no API key
        try:
            url = f"{self.base_url}/records/{db_id}/details?fields=Formula,Mol3d" # Request Mol3d
            response = self.session.request("GET", url, headers=self.headers, timeout=30)
            response.raise_for_status()
            data = response.json()
            mol3d_content = data.get('mol3d')

            if mol3d_content:
                # Basic check if it looks like MOL/SDF format (very rudimentary)
                 if "V2000" in mol3d_content or "V3000" in mol3d_content:
                     # Convert string content to bytes
                     return mol3d_content.encode('utf-8')
                 else:
                     logging.warning(f"ChemSpider download for {name} (ID: {db_id}) did not return expected Mol3d format.")
                     return None
            else:
                 logging.debug(f"No Mol3d content found for {name} (ID: {db_id}) in ChemSpider details.")
                 return None

        except requests.exceptions.RequestException as e:
            status_code = e.response.status_code if e.response is not None else 'N/A'
            logging.warning(f"ChemSpider download HTTP error for {name} (ID: {db_id}, Status: {status_code}): {e}")
            return None
        except json.JSONDecodeError as e:
            logging.error(f"ChemSpider download JSON decode error for {name} (ID: {db_id}): {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected ChemSpider download error for {name} (ID: {db_id}): {e}")
            return None

class MultiDatabaseDownloader:
    """Downloads 3D structures from multiple chemical databases."""
    def __init__(self, chemspider_api_key: Optional[str], max_workers: int = 4):
        self.max_workers = max_workers
        self.chemspider_api_key = chemspider_api_key
        self._setup_logging()
        self._initialize_databases()
        self.downloaded_files: Dict[str, Path] = {} # Store name -> filepath mapping
        self.failed_compounds: Set[str] = set()

    def _setup_logging(self):
        log_dir = Path("logs")
        log_dir.mkdir(exist_ok=True)
        self.log_file = log_dir / 'structure_download.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s [%(levelname)s] (%(threadName)s) %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler() # Also print to console
            ]
        )
        logging.info("--- Starting new download session ---")

    def _initialize_databases(self):
        self.databases = {
            'pubchem': PubChemDB(),
            # Only add ChemSpider if API key is provided
            'chemspider': ChemSpiderDB(self.chemspider_api_key) if self.chemspider_api_key else None,
        }
        # Filter out None databases
        self.databases = {k: v for k, v in self.databases.items() if v is not None}
        if not self.databases:
             raise ValueError("No chemical databases could be initialized. Check API keys and configurations.")
        logging.info(f"Initialized databases: {list(self.databases.keys())}")

    def _clean_filename(self, name: str) -> str:
        """Creates a filesystem-safe filename from a compound name."""
        safe_name = str(name)
        # Replace common problematic characters
        invalid_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '|', '+', '%', 'and']
        for char in invalid_chars:
            safe_name = safe_name.replace(char, '_')
        # Replace multiple spaces/underscores with a single underscore
        safe_name = re.sub(r'[\s_]+', '_', safe_name)
        # Remove leading/trailing spaces, dots, underscores
        safe_name = safe_name.strip('. _')
        # Limit length (e.g., 100 chars) to avoid issues on some systems
        safe_name = safe_name[:100]
        return safe_name if safe_name else "unknown_compound"

    def save_structure(self, content: bytes, search_name: str, output_dir: Path) -> Optional[Path]:
        """Saves structure content to a unique file in the specified directory."""
        try:
            safe_name = self._clean_filename(search_name)
            base_filepath = output_dir / f"{safe_name}.sdf"
            output_filepath = base_filepath
            counter = 1
            # Ensure unique filename if it already exists
            while output_filepath.exists():
                output_filepath = output_dir / f"{safe_name}_{counter}.sdf"
                counter += 1

            output_dir.mkdir(parents=True, exist_ok=True) # Ensure directory exists
            with open(output_filepath, 'wb') as f:
                f.write(content)
            logging.info(f"Saved structure for '{search_name}' to {output_filepath}")
            return output_filepath
        except IOError as e:
            logging.error(f"I/O error saving structure for '{search_name}': {e}")
            return None
        except Exception as e:
            logging.error(f"Unexpected error saving structure for '{search_name}': {e}")
            return None

    def download_structure(self, name: str, output_dir: Path) -> Optional[Path]:
        """Tries to download structure from configured databases."""
        logging.debug(f"Attempting to download structure for: '{name}'")
        # Prioritize PubChem as it often has better 3D conformers
        db_priority = ['pubchem', 'chemspider']

        for db_name in db_priority:
            database = self.databases.get(db_name)
            if not database:
                continue # Skip if this database wasn't initialized

            logging.info(f"Trying {db_name.capitalize()} for '{name}'")
            db_id = database.search_compound(name)
            if db_id:
                logging.info(f"Found {db_name.capitalize()} ID {db_id} for '{name}'")
                content = database.download_3d_structure(db_id, name)
                if content:
                    saved_path = self.save_structure(content, name, output_dir)
                    if saved_path:
                        logging.info(f"Successfully downloaded and saved '{name}' from {db_name.capitalize()}.")
                        return saved_path # Success! Return the path
                    else:
                         logging.error(f"Failed to save structure file for '{name}' from {db_name.capitalize()}.")
                         # Continue to next database if saving failed
                else:
                    logging.warning(f"Found ID but failed to download 3D structure for '{name}' from {db_name.capitalize()}.")
            else:
                 logging.info(f"Compound '{name}' not found in {db_name.capitalize()}.")

        logging.warning(f"Failed to download '{name}' from any configured database.")
        self.failed_compounds.add(name)
        return None

    def process_compound_list(self, compounds: List[str], output_dir_name: str):
        """Processes a list of compounds, downloading structures to a specific directory."""
        specific_output_dir = Path(output_dir_name)
        specific_output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"--- Processing list for output directory: {specific_output_dir} ---")

        self.downloaded_files.clear()
        self.failed_compounds.clear()
        processed_count = 0

        with ThreadPoolExecutor(max_workers=self.max_workers, thread_name_prefix='Downloader') as executor:
            # Map compound name to future object
            future_to_compound = {
                executor.submit(self.download_structure, compound, specific_output_dir): compound
                for compound in set(compounds) # Process unique compounds
            }

            total_compounds = len(future_to_compound)
            with tqdm(total=total_compounds, desc=f"Downloading for {output_dir_name}", unit="compound") as pbar:
                for future in as_completed(future_to_compound):
                    compound_name = future_to_compound[future]
                    processed_count += 1
                    try:
                        result_path = future.result() # Returns Path if successful, None otherwise
                        if result_path:
                            self.downloaded_files[compound_name] = result_path
                        # Failure is handled within download_structure by adding to self.failed_compounds
                    except Exception as e:
                        logging.error(f"Thread error processing '{compound_name}': {e}")
                        self.failed_compounds.add(compound_name) # Mark as failed on thread error
                    finally:
                         pbar.update(1) # Update progress bar regardless of outcome

        logging.info(f"--- Finished processing for {specific_output_dir} ---")
        successfully_downloaded = len(self.downloaded_files)
        failed_count = len(self.failed_compounds)
        total_processed_unique = successfully_downloaded + failed_count

        logging.info(f"Total unique compounds attempted: {total_processed_unique}")
        logging.info(f"Successfully downloaded: {successfully_downloaded}")
        logging.info(f"Failed to download: {failed_count}")

        failed_list_filename = specific_output_dir.parent / f"{specific_output_dir.name}_failed_downloads.txt"
        try:
            with open(failed_list_filename, 'w', encoding='utf-8') as f:
                if self.failed_compounds:
                    f.write(f"# Compounds that could not be downloaded for list processed into '{specific_output_dir.name}'\n")
                    for compound in sorted(list(self.failed_compounds)):
                        f.write(f"{compound}\n")
                else:
                     f.write("# No compounds failed to download for this list.\n")
            logging.info(f"List of failed compounds saved to: {failed_list_filename}")
        except IOError as e:
             logging.error(f"Could not write failed compounds list {failed_list_filename}: {e}")

        return {
            "total_attempted": total_processed_unique,
            "downloaded": successfully_downloaded,
            "failed": failed_count,
            "output_dir": str(specific_output_dir),
            "failed_list_file": str(failed_list_filename)
        }
def main():
    """Configures and runs the structure download process for multiple compound lists."""
    # --- Configuration ---
    try:
        with open('config.json', 'r') as f:
            config_data = json.load(f)
            # Check specifically for ChemSpider key
            chemspider_api_key = config_data.get('CHEMSPIDER_API_KEY')
            if not chemspider_api_key:
                 print("Warning: 'CHEMSPIDER_API_KEY' not found in config.json. ChemSpider searches will be skipped.")
    except FileNotFoundError:
        print("Error: config.json file not found.")
        print("Please create config.json. It should contain at least:")
        print('{\n  "CHEMSPIDER_API_KEY": "YOUR_API_KEY_HERE_OR_NULL"\n}')
        return
    except json.JSONDecodeError:
        print("Error: config.json is not valid JSON.")
        return

    compound_lists_config = [
        {"input_file": "lippia_phytochemicals_unique_list.txt", "output_dir": "3d_structures_lippia"},
        {"input_file": "Atermisia_phytochemicals_unique_list.txt", "output_dir": "3d_structures_artemisia"},
        {"input_file": "allium_phytochemicals_unique_list.txt", "output_dir": "3d_structures_allium"},
    ]

    base_output_dir = Path("downloaded_structures") # Parent directory for all results
    base_output_dir.mkdir(exist_ok=True)
    try:
         downloader = MultiDatabaseDownloader(
             chemspider_api_key=chemspider_api_key,
             max_workers=4 # Adjust number of workers as needed
         )
    except ValueError as e: # Catches initialization errors (e.g., no usable databases)
        print(f"Error initializing downloader: {e}")
        return

    all_run_stats = [] # Store stats from each run

    for config in compound_lists_config:
        input_file = Path(config["input_file"])
        output_dir = base_output_dir / config["output_dir"]

        print(f"\n===== Processing: {input_file} =====")
        if not input_file.exists():
            print(f"Error: Input file not found: {input_file}")
            logging.error(f"Input file not found: {input_file}, skipping.")
            all_run_stats.append({**config, "status": "skipped - input file not found"})
            continue

        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                # Read lines, strip whitespace, filter out empty lines
                compounds = [line.strip() for line in f if line.strip()]
                if not compounds:
                     print(f"Warning: Input file {input_file} is empty or contains only whitespace.")
                     logging.warning(f"Input file {input_file} is empty, skipping.")
                     all_run_stats.append({**config, "status": "skipped - input file empty"})
                     continue

            print(f"Found {len(set(compounds))} unique compounds in {input_file}.")

            run_stats = downloader.process_compound_list(compounds, str(output_dir))
            all_run_stats.append({**config, **run_stats, "status": "completed"})

        except IOError as e:
            print(f"Error reading input file {input_file}: {e}")
            logging.error(f"Error reading input file {input_file}: {e}, skipping.")
            all_run_stats.append({**config, "status": f"skipped - error reading input file: {e}"})
        except Exception as e:
            print(f"An unexpected error occurred while processing {input_file}: {e}")
            logging.exception(f"Unexpected error processing {input_file}, skipping.") # Log full traceback
            all_run_stats.append({**config, "status": f"skipped - unexpected error: {e}"})

    print("\n\n===== Overall Download Summary =====")
    total_attempted_all = 0
    total_downloaded_all = 0
    total_failed_all = 0

    for stats in all_run_stats:
        print(f"\n--- Summary for: {stats['input_file']} ---")
        print(f"  Status: {stats['status']}")
        if stats['status'] == 'completed':
            print(f"  Output Directory: {stats['output_dir']}")
            print(f"  Compounds Attempted: {stats['total_attempted']}")
            print(f"  Successfully Downloaded: {stats['downloaded']}")
            print(f"  Failed to Download: {stats['failed']}")
            print(f"  Failed List File: {stats['failed_list_file']}")
            total_attempted_all += stats['total_attempted']
            total_downloaded_all += stats['downloaded']
            total_failed_all += stats['failed']
        else:
            print(f"  Processing was skipped.")

    print("\n--- Grand Totals ---")
    print(f"Total unique compounds attempted across all lists: {total_attempted_all}")
    print(f"Total structures successfully downloaded across all lists: {total_downloaded_all}")
    print(f"Total unique compounds failed across all lists: {total_failed_all}")
    print(f"\nCheck '{downloader.log_file}' for detailed logging information.")
    print("===================================")

if __name__ == "__main__":
    import re
    main()
