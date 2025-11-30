import argparse
import csv
import hashlib
import yaml
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter, Retry


def session_with_retry():
    retries = Retry(
        total=5,
        backoff_factor=0.5,
        status_forcelist=[429, 500, 502, 503, 504],
    )
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retries))
    return s


def verify_md5(path, md5):
    if not md5:
        return True
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest() == md5


def download_file(url, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    session = session_with_retry()

    if out_path.exists():
        print(f"[SKIP] Exists: {out_path}")
        return

    print(f"[DOWNLOAD] {url} → {out_path}")
    with session.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1_048_576):
                if chunk:
                    f.write(chunk)


def run_download(config):
    base = Path(config["base_dir"])
    raw_dir = base / config["raw_data_dir"]
    csv_path = Path("scripts/download_links.csv")

    print(f"[INFO] Reading download list: {csv_path}")
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            file_id = row["id"]
            url = row["url"]
            md5 = row.get("md5", None)

            out_path = raw_dir / f"{file_id}"

            download_file(url, out_path)

            if md5 and not verify_md5(out_path, md5):
                print(f"[ERROR] MD5 mismatch: {file_id}")
            else:
                print(f"[OK] {file_id} downloaded.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    run_download(config)


if __name__ == "__main__":
    main()
