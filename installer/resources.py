from pathlib import Path
from .utils import ensure_dir, download_file, log


def install_genome_resources(genome: str, genome_entry: dict, resource_root: str, force_redownload: bool = False):
    genome_resource_dir = Path(resource_root) / genome
    ensure_dir(genome_resource_dir)

    resource_payload = {}
    download_records = []

    resources = genome_entry.get("resources", {})
    for name, info in resources.items():
        url = info["url"]
        filename = info["filename"]
        output_path = genome_resource_dir / filename
        saved_path = download_file(url, str(output_path), force=force_redownload)
        resource_payload[name] = saved_path
        download_records.append({
            "name": name,
            "genome": genome,
            "url": url,
            "path": saved_path
        })
        log(f"Installed resource for {genome}: {name} -> {saved_path}")

    return str(genome_resource_dir.resolve()), resource_payload, download_records
