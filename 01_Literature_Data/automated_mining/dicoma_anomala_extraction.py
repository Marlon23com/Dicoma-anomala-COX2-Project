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


class PubMedDownloader:
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def __init__(self,
                 search_term: str,
                 download_dir_suffix: str,
                 queries: Optional[List[str]] = None,
                 require_review: bool = False,
                 require_phytochem: bool = False,
                 include_extra_metadata: bool = False,
                 max_results_per_query: int = 100,
                 max_filename_length: int = 50,
                 base_download_dir: str = '.',
                 min_year: Optional[int] = None,
                 max_year: Optional[int] = None):

        self.search_term = search_term
        self.download_dir = os.path.join(
            base_download_dir,
            self.clean_filename(f"{search_term}_{download_dir_suffix}", max_length=100)
        )
        self.queries = queries if queries else [self.search_term]
        self.require_review = require_review
        self.require_phytochem = require_phytochem
        self.include_extra_metadata = include_extra_metadata
        self.max_results_per_query = max_results_per_query
        self.max_filename_length = max_filename_length
        self.min_year = min_year
        self.max_year = max_year

        os.makedirs(self.download_dir, exist_ok=True)
        print(f"Initialized downloader for '{search_term}'. Saving to: {self.download_dir}")
        print(f"  Require Review: {self.require_review}, Require Phytochem: {self.require_phytochem}, Extra Metadata: {self.include_extra_metadata}")

    def clean_filename(self, text: str, max_length: int) -> str:
        if not text:
            text = 'Untitled'
        clean_text = re.sub(r'[^\w\s-]', '', text)
        clean_text = re.sub(r'\s+', '_', clean_text)
        clean_text = clean_text.strip('_-')
        return clean_text[:max_length]

    def _fetch_pmids(self, query: str) -> List[str]:
        search_query_parts = [f"({query})"]

        if self.require_review:
            search_query_parts.append("(review[Publication Type])")

        if self.require_phytochem:
            search_query_parts.append("(phytochemic* OR chemical constituent* OR compound* OR metabolite*)")

        if self.min_year is not None and self.max_year is not None:
            search_query_parts.append(f"(\"{self.min_year}\"[Date - Publication] : \"{self.max_year}\"[Date - Publication])")

        full_query = " AND ".join(search_query_parts)
        search_url = f"{self.BASE_URL}esearch.fcgi?db=pubmed&term={quote_plus(full_query)}&retmax={self.max_results_per_query}&retmode=xml"

        print(f"  Searching PubMed with query: {full_query}")
        try:
            response = requests.get(search_url, timeout=30)
            response.raise_for_status()
            root = ET.fromstring(response.content)
            pmids = [id_elem.text for id_elem in root.findall('.//IdList/Id') if id_elem.text]
            print(f"  Found {len(pmids)} PMIDs for query.")
            return pmids
        except requests.exceptions.RequestException as e:
            print(f"  PubMed search network error for query '{query}': {e}")
            return []
        except ET.ParseError as e:
            print(f"  PubMed search XML parsing error for query '{query}': {e}")
            print(f"  Response content: {response.text[:500]}...")
            return []
        except Exception as e:
            print(f"  Unexpected PubMed search error for query '{query}': {e}")
            return []

    def _fetch_article_details(self, pmid: str) -> Optional[Dict]:
        fetch_url = f"{self.BASE_URL}efetch.fcgi?db=pubmed&id={pmid}&retmode=xml"
        try:
            time.sleep(random.uniform(0.4, 1.0))
            fetch_response = requests.get(fetch_url, timeout=30)
            fetch_response.raise_for_status()
            article_root = ET.fromstring(fetch_response.content)

            title_elem = article_root.find('.//ArticleTitle')
            abstract_parts = article_root.findall('.//Abstract/AbstractText')
            full_abstract = "\n".join([part.text for part in abstract_parts if part.text]) if abstract_parts else None

            article_info = {
                'pmid': pmid,
                'title': title_elem.text.strip() if title_elem is not None and title_elem.text else 'No Title Found',
                'abstract': full_abstract if full_abstract else 'No Abstract Available',
                'publication_types': [elem.text for elem in article_root.findall('.//PublicationType') if elem.text]
            }

            if self.require_review and not any('Review' in pt for pt in article_info['publication_types']):
                return None

            if self.include_extra_metadata:
                year_elem = article_root.find('.//PubDate/Year')
                article_info['year'] = year_elem.text if year_elem is not None else 'Year Unknown'

                author_list = []
                for author in article_root.findall('.//Author'):
                    last_name = author.find('LastName')
                    fore_name = author.find('ForeName')
                    initials = author.find('Initials')
                    if last_name is not None and last_name.text:
                        name = f"{last_name.text}, {fore_name.text}" if fore_name is not None and fore_name.text else f"{last_name.text}, {initials.text}" if initials is not None and initials.text else last_name.text
                        author_list.append(name)
                article_info['authors'] = '; '.join(author_list) if author_list else 'Authors Unknown'

                keywords = set()
                for k_elem in article_root.findall('.//KeywordList/Keyword'):
                    if k_elem.text: keywords.add(k_elem.text.strip())
                for mesh_elem in article_root.findall('.//MeshHeadingList/MeshHeading/DescriptorName'):
                    if mesh_elem.text: keywords.add(mesh_elem.text.strip())
                article_info['keywords'] = '; '.join(sorted(list(keywords))) if keywords else 'No Keywords Found'

            return article_info

        except requests.exceptions.RequestException as e:
            print(f"    Network error fetching details for PMID {pmid}: {e}")
            return None
        except ET.ParseError as e:
            print(f"    XML parsing error for PMID {pmid}: {e}")
            print(f"    Response content: {fetch_response.text[:500]}...")
            return None
        except Exception as e:
            print(f"    Unexpected error fetching details for PMID {pmid}: {e}")
            return None

    def download_article_content(self, article: Dict) -> Optional[str]:
        try:
            safe_title = self.clean_filename(article.get('title', 'Untitled'), self.max_filename_length)
            pmid = article.get('pmid', 'unknown_pmid')
            filename = os.path.join(self.download_dir, f"{safe_title}_{pmid}.txt")

            with open(filename, 'w', encoding='utf-8') as f:
                f.write(f"Title: {article.get('title', 'N/A')}\n\n")
                if 'authors' in article:
                    f.write(f"Authors: {article.get('authors', 'N/A')}\n\n")
                if 'year' in article:
                    f.write(f"Year: {article.get('year', 'N/A')}\n\n")
                f.write(f"PMID: {pmid}\n\n")
                if 'keywords' in article:
                    f.write(f"Keywords/MeSH: {article.get('keywords', 'N/A')}\n\n")
                if article.get('publication_types'):
                    f.write(f"Publication Types: {', '.join(article['publication_types'])}\n\n")
                f.write(f"Abstract:\n{article.get('abstract', 'N/A')}\n")

            return filename
        except IOError as e:
            print(f"  Error writing file for PMID {article.get('pmid', 'N/A')}: {e}")
            return None
        except Exception as e:
            print(f"  Unexpected error saving article for PMID {article.get('pmid', 'N/A')}: {e}")
            return None

    def run(self) -> List[str]:
        print(f"\n--- Starting download process for: {self.search_term} ---")
        downloaded_files: List[str] = []
        processed_pmids: Set[str] = set()

        for query in self.queries:
            print(f"\nProcessing query: '{query}'")
            pmids = self._fetch_pmids(query)

            if not pmids:
                print("  No PMIDs found or error occurred for this query.")
                continue

            count = 0
            for pmid in pmids:
                if pmid in processed_pmids:
                    continue

                print(f"  Fetching details for PMID: {pmid}...")
                article_details = self._fetch_article_details(pmid)
                processed_pmids.add(pmid)

                if article_details:
                    saved_path = self.download_article_content(article_details)
                    if saved_path:
                        print(f"    Successfully saved: {os.path.basename(saved_path)}")
                        downloaded_files.append(saved_path)
                        count += 1
                    else:
                        print(f"    Failed to save article content for PMID {pmid}.")

            print(f"  Downloaded {count} new articles for this query.")
            time.sleep(random.uniform(1, 3))

        total_downloaded = len(downloaded_files)
        print(f"\n--- Download process finished for: {self.search_term} ---")
        print(f"Total unique articles downloaded: {total_downloaded}")
        return downloaded_files


class UnifiedScopusSearcher:
    def __init__(self,
                 species_name: str,
                 phyto_keywords: List[str],
                 download_dir: str,
                 api_key: str):

        self.species_name = species_name
        self.phyto_keywords = phyto_keywords
        self.download_dir = download_dir
        self.api_key = api_key

        if not api_key:
            raise ValueError("Scopus API key is required")

        os.makedirs(self.download_dir, exist_ok=True)
        self.setup_client()

        self.doc_type_map = {
            'ar': 'Article', 're': 'Review', 'cp': 'Conference Paper',
            'ch': 'Book Chapter', 'bk': 'Book', 'le': 'Letter',
            'no': 'Note', 'sh': 'Short Survey'
        }
        self.all_results = []

    def setup_client(self):
        try:
            self.els_client = ElsClient(self.api_key)
            test_search = ElsSearch('TITLE(test)', 'scopus')
            test_search.execute(self.els_client, get_all=False)
            print(f"Scopus client configured successfully for {self.species_name} - API key is valid.")
        except Exception as e:
            print(f"Error configuring Scopus client for {self.species_name}: {e}")
            raise

    def _build_query(self, min_year: int) -> str:
        query_parts = [
            f'TITLE-ABS-KEY("{self.species_name}")',
            f'PUBYEAR > {min_year}',
            'LANGUAGE(english)'
        ]

        if self.phyto_keywords:
            keyword_query_part = ' OR '.join(f'"{k}"' for k in self.phyto_keywords)
            combined_phyto_query = (
                f'(TITLE-ABS-KEY(phytochemical OR "chemical composition") OR '
                f'TITLE-ABS-KEY({keyword_query_part}))'
            )
            query_parts.append(combined_phyto_query)
        else:
            query_parts.append('TITLE-ABS-KEY(phytochemical OR "chemical composition")')

        return " AND ".join(query_parts)

    def search_scopus(self, max_results: int = 2000, min_year: int = 2000) -> List[Dict]:
        query = self._build_query(min_year)
        print(f"\nSearching Scopus for: {self.species_name}")
        print(f"Executing search with query: {query}")

        try:
            search = ElsSearch(query, 'scopus')
            search.execute(self.els_client, get_all=True)

            publications = []
            results_to_process = search.results[:max_results]
            print(f"Fetched {len(search.results)} total results, processing up to {max_results}.")

            for result in results_to_process:
                doc_type_code = result.get('prism:aggregationType', '').lower()
                doc_type = self.doc_type_map.get(doc_type_code, doc_type_code.capitalize())

                pub_data = {
                    'title': result.get('dc:title', 'N/A'),
                    'abstract': result.get('dc:description', 'N/A'),
                    'authors': [author.get('authname', '') for author in result.get('author', []) if author.get('authname')],
                    'journal': result.get('prism:publicationName', 'N/A'),
                    'year': result.get('prism:coverDate', 'N/A')[:4],
                    'doi': result.get('prism:doi', 'N/A'),
                    'cited_by': result.get('citedby-count', '0'),
                    'scopus_id': result.get('dc:identifier', '').replace('SCOPUS_ID:', ''),
                    'document_type': doc_type,
                    'keywords': [kw.strip() for kw in result.get('authkeywords', '').split('|') if kw.strip()] if result.get('authkeywords') else [],
                    'source': 'Scopus'
                }
                publications.append(pub_data)

            print(f"Successfully processed {len(publications)} publications for {self.species_name}")
            self.all_results = publications
            return publications

        except Exception as e:
            print(f"Detailed error during Scopus search for {self.species_name}:")
            print(f"Error message: {str(e)}")
            traceback.print_exc()
            return []

    def save_results(self) -> Optional[str]:
        if not self.all_results:
            print(f"No results to save for {self.species_name}.")
            return None

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        safe_species_name = self.species_name.replace(' ', '_').lower()
        base_filename = f"{safe_species_name}_{timestamp}"

        json_filename = os.path.join(self.download_dir, f"{base_filename}_detailed.json")
        try:
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(self.all_results, f, indent=2, ensure_ascii=False)
        except IOError as e:
            print(f"Error saving JSON file {json_filename}: {e}")
            return None

        excel_filename = os.path.join(self.download_dir, f"{base_filename}_summary.xlsx")
        try:
            summary_data = [{
                'Title': pub.get('title', 'N/A'),
                'Year': pub.get('year', 'N/A'),
                'Document Type': pub.get('document_type', 'N/A'),
                'Journal': pub.get('journal', 'N/A'),
                'Authors': '; '.join(pub.get('authors', [])),
                'Keywords': '; '.join(pub.get('keywords', [])),
                'DOI': pub.get('doi', 'N/A'),
                'Citations': int(pub.get('cited_by', 0)),
                'Scopus ID': pub.get('scopus_id', 'N/A')
            } for pub in self.all_results]

            df = pd.DataFrame(summary_data)
            doc_type_summary = df['Document Type'].value_counts()

            with pd.ExcelWriter(excel_filename) as writer:
                df.to_excel(writer, sheet_name='Publications', index=False)
                doc_type_summary.to_frame('Count').to_excel(writer, sheet_name='Document Type Summary')

        except (IOError, ImportError) as e:
            print(f"Error saving Excel file {excel_filename}: {e}")
            print("Ensure 'pandas' and 'openpyxl' are installed.")
        except Exception as e:
            print(f"Unexpected error saving Excel file {excel_filename}: {e}")

        print(f"\nResults for {self.species_name} saved to:")
        print(f"- JSON: {json_filename}")
        if os.path.exists(excel_filename):
            print(f"- Excel: {excel_filename}")

        return base_filename

    def run(self, max_results: int = 2000, min_year: int = 2000) -> Optional[str]:
        self.search_scopus(max_results=max_results, min_year=min_year)
        base_filename = self.save_results()

        print(f"\n--- Results Summary for {self.species_name} ---")
        print(f"Total results processed: {len(self.all_results)}")
        doc_type_counts = defaultdict(int)
        for result in self.all_results:
            doc_type_counts[result.get('document_type', 'Unknown')] += 1
        if doc_type_counts:
            print("\nDocument type breakdown:")
            for doc_type, count in sorted(doc_type_counts.items()):
                print(f"- {doc_type}: {count}")
        print("--- End Summary ---")

        return base_filename


def main():
    SEARCH_QUERIES = [
        # 1. Core plant + phytochemistry
        '("Dicoma anomala" OR "Chifumuro" OR "Dicoma") AND ("phytochemical*" OR "constituent*" OR "metabolite*" OR "isolated compound*" OR "chemical composition")',
        
        # 2. Plant + traditional/respiratory use
        '("Dicoma anomala" OR "Chifumuro") AND ("respiratory" OR "cough" OR "cold" OR "fever" OR "nasal congestion" OR "anti-inflammatory" OR "antimicrobial")',
        
        # 3. Target compounds (specific isolates)
        '("Dehydrobrachylaenolide" OR "Dicomanol" OR "Lupeol" OR "Lupenone" OR "Beta-amyrin" OR "Quercetin" OR "Luteolin") AND ("Dicoma" OR "Asteraceae" OR "anti-inflammatory")',
        
        # 4. COX-2/Sinusitis context + natural products (for background & justification)
        '("Chronic rhinosinusitis" OR "sinusitis" OR "nasal inflammation") AND ("COX-2" OR "cyclooxygenase-2" OR "PTGS2") AND ("natural product*" OR "plant extract*" OR "phytochemical*")'
    ]

    pubmed_downloader = PubMedDownloader(
        search_term="Dicoma anomala",
        download_dir_suffix="dicoma_anomala_articles",
        queries=SEARCH_QUERIES,
        require_review=False,
        require_phytochem=True,
        include_extra_metadata=True,
        max_results_per_query=200,
        min_year=2000,
        max_year=2026
    )
    pubmed_files = pubmed_downloader.run()

    try:
        with open('api.txt', 'r') as file:
            api_key = file.read().strip()
        if not api_key:
            raise ValueError("API key in api.txt is empty.")
    except FileNotFoundError:
        print("Error: api.txt file not found in the current directory.")
        print("Please create api.txt and paste your Scopus API key into it.")
        api_key = None
    except ValueError as e:
        print(f"Error: {e}")
        api_key = None

    all_base_filenames = {}

    if api_key:
        base_output_dir = 'scopus_dicoma_anomala_results'
        os.makedirs(base_output_dir, exist_ok=True)

        try:
            searcher = UnifiedScopusSearcher(
                species_name="Dicoma anomala",
                phyto_keywords=[
                    'sesquiterpene', 'sesquiterpene lactone', 'flavonoid',
                    'phenolic compound', 'anti-inflammatory', 'COX-2',
                    'chemical constituent', 'phytochemical', 'essential oil'
                ],
                download_dir=os.path.join(base_output_dir, 'dicoma_anomala_publications'),
                api_key=api_key
            )
            base_filename = searcher.run(max_results=2000, min_year=2000)
            if base_filename:
                all_base_filenames["Dicoma anomala"] = base_filename
        except Exception as e:
            print(f"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print(f"A critical error occurred processing Dicoma anomala: {e}")
            print("Check API key validity and network connection.")
            traceback.print_exc()
            print(f"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    print("\n\n===== Dicoma anomala Extraction Summary =====")
    print(f"PubMed articles downloaded: {len(pubmed_files)}")
    if all_base_filenames:
        print("Scopus results saved with the following base filenames:")
        for species, fname in all_base_filenames.items():
            print(f"- {species}: {fname}")
    else:
        if api_key:
            print("No Scopus results were successfully saved.")
        else:
            print("Scopus search skipped because api.txt was missing or invalid.")
    print("===========================================")


if __name__ == "__main__":
    main()
