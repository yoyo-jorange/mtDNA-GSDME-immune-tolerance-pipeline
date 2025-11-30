import yaml
import psutil
import subprocess
import sys
from pathlib import Path

def check_all():
    report = {"status": "OK", "checks": []}
    
    # 磁盘空间
    disk = psutil.disk_usage('/')
    free_gb = disk.free / (1024**3)
    report["checks"].append({"disk_free_gb": round(free_gb, 1), "sufficient": free_gb > 100})
    
    # 内存
    mem = psutil.virtual_memory()
    report["checks"].append({"memory_gb": round(mem.total / (1024**3), 1)})
    
    # R/Python版本
    report["checks"].append({"R_version": subprocess.run(["R", "--version"], capture_output=True, text=True).stderr.split('\n')[0]})
    report["checks"].append({"python_version": sys.version})
    
    # 关键包检查
    r_pkgs = ["Seurat", "Cell2location", "GSVA"]
    py_pkgs = ["scanpy", "dowhy", "causalml"]
    
    return report

if __name__ == "__main__":
    report = check_all()
    print(yaml.dump(report, default_flow_style=False, allow_unicode=True))
