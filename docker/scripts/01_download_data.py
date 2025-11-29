import os
import csv
import time
import subprocess
from pathlib import Path

RAW_DIR = Path("data_raw")
LINKS_CSV = Path("scripts/download_links.csv")

def download_file(url, out_path, max_retries=3):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, max_retries + 1):
        print(f"下载 {url} → {out_path} (第 {attempt} 次尝试)")
        try:
            # 使用 curl 下载
            result = subprocess.run(
                ["curl", "-L", "-o", str(out_path), url],
                check=False,
            )
            if result.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0:
                print(f"完成: {out_path}")
                return True
        except Exception as e:
            print(f"错误: {e}")
        time.sleep(3)
    print(f"[失败] 无法下载: {url}")
    return False

def main():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    log_path = Path("results/tables/download_summary.csv")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    results = []
    if not LINKS_CSV.exists():
        print(f"[警告] 找不到 {LINKS_CSV}，请让学生/助手填写下载链接表格。")
        return

    with open(LINKS_CSV, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            url = row["url"]
            rel_path = row["rel_path"]  # 例如 "10x/visium_hd/sample1.tar.gz"
            out_path = RAW_DIR / rel_path
            ok = download_file(url, out_path)
            results.append({"url": url, "rel_path": rel_path, "status": "ok" if ok else "fail"})

    with open(log_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["url", "rel_path", "status"])
        writer.writeheader()
        writer.writerows(results)

    print(f"下载总结已保存: {log_path}")

if __name__ == "__main__":
    main()
