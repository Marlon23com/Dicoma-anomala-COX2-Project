import os
import requests
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import time

class SDFDownloader:
    def __init__(self, base_dir="02_Ligands/raw_sdf", max_workers=5, delay=0.1):
        self.base_dir = base_dir
        self.max_workers = max_workers
        self.delay = delay
        self.session = requests.Session()
        # Set user agent to avoid blocking
        self.session.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'})

        # Ensure base directory exists
        os.makedirs(self.base_dir, exist_ok=True)

        # Setup logging
        logging.basicConfig(filename='logs/sdf_download.log', level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')

    def download_sdf(self, cid, compound_name, subfolder):
        """Download 3D SDF for a given CID"""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type=3d"
        output_dir = os.path.join(self.base_dir, subfolder)
        os.makedirs(output_dir, exist_ok=True)

        safe_name = "".join(c for c in compound_name if c.isalnum() or c in (' ', '-', '_')).rstrip()
        filename = f"{safe_name}.sdf"
        filepath = os.path.join(output_dir, filename)

        try:
            response = self.session.get(url, timeout=30)
            if response.status_code == 200:
                with open(filepath, 'w') as f:
                    f.write(response.text)
                logging.info(f"Downloaded {filename} for CID {cid}")
                return True
            else:
                logging.warning(f"Failed to download CID {cid}: HTTP {response.status_code}")
                return False
        except Exception as e:
            logging.error(f"Error downloading CID {cid}: {str(e)}")
            return False
        finally:
            time.sleep(self.delay)  # Rate limiting

    def process_cid_list(self, excel_file, subfolder_name):
        """Process CIDs from an Excel file"""
        try:
            df = pd.read_excel(excel_file)
            if 'cid' not in df.columns:
                logging.error(f"No 'cid' column in {excel_file}")
                return

            cids_and_names = []
            for _, row in df.iterrows():
                cid = row['cid']
                name = row['Compound Name']
                if pd.notna(cid) and pd.notna(name):
                    cids_and_names.append((int(cid), str(name)))

            logging.info(f"Found {len(cids_and_names)} compounds in {excel_file}")

            failed_downloads = []
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = [executor.submit(self.download_sdf, cid, name, subfolder_name)
                          for cid, name in cids_and_names]

                for future in tqdm(as_completed(futures), total=len(futures), desc=f"Downloading {subfolder_name}"):
                    if not future.result():
                        # Note: we don't know which one failed, but we can log it later if needed
                        pass

            # Save failed downloads if any
            if failed_downloads:
                failed_file = os.path.join(self.base_dir, subfolder_name, "failed_downloads.txt")
                with open(failed_file, 'w') as f:
                    for cid in failed_downloads:
                        f.write(f"{cid}\n")

        except Exception as e:
            logging.error(f"Error processing {excel_file}: {str(e)}")

if __name__ == "__main__":
    downloader = SDFDownloader()

    excel_files = [
        ("01_Literature_Data/manual_curation/endophyte_derived_for_d_anomala.csv.xlsx", "endophyte_derived_for_d_anomala"),
        ("01_Literature_Data/manual_curation/high_confidence_compounds.csv.xlsx", "high_confidence_compounds"),
        ("01_Literature_Data/manual_curation/require_structural_curation.csv.xlsx", "require_structural_curation")
    ]

    for excel_file, subfolder in excel_files:
        if os.path.exists(excel_file):
            downloader.process_cid_list(excel_file, subfolder)
        else:
            logging.error(f"Excel file not found: {excel_file}")