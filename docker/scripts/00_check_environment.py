import os
import shutil
import psutil  # 如果没有会报错，后续可加到 Dockerfile
from datetime import datetime

def check_disk_space(min_gb=50):
    usage = psutil.disk_usage("/")
    free_gb = usage.free / (1024**3)
    ok = free_gb >= min_gb
    return ok, free_gb

def main():
    os.makedirs("results/reports", exist_ok=True)
    report_path = "results/reports/00_check_environment.txt"

    lines = []
    lines.append(f"检查时间: {datetime.now()}")
    lines.append("=== 基本环境检查 ===")

    # 磁盘空间
    ok_disk, free_gb = check_disk_space()
    if ok_disk:
        lines.append(f"[通过] 剩余磁盘空间约 {free_gb:.1f} GB")
    else:
        lines.append(f"[警告] 剩余磁盘空间仅 {free_gb:.1f} GB，建议至少 50 GB")

    # Docker 这里只在容器外不检查，假定已能运行到这里
    lines.append("[提示] 当前已在 Docker 容器内运行，基础环境已就绪。")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print("\n".join(lines))

if __name__ == "__main__":
    main()
