#!/bin/bash

#===============================================================================
# SCRIPT: convert_homer_html_to_pdf.sh
# PURPOSE: Convert HTML reports to PDF format
#
# DESCRIPTION:
# Converts HTML reports to PDF format for easier viewing and sharing.
#
# REQUIREMENTS:
# - Python with weasyprint package (installed automatically if needed)
#
# USAGE:
# ./scripts/convert_homer_html_to_pdf.sh              # Convert all default files
# sbatch scripts/convert_homer_html_to_pdf.sh         # Submit to SLURM
#
# INPUTS:
# - results/multiqc/summary_report.html
# - results/multiqc/TES_TEAD1_CutTag_Report.html
# - results/report/TES_TEAD1_CutTag_Report.html
#
# OUTPUTS:
# - results/multiqc/summary_report.pdf
# - results/multiqc/TES_TEAD1_CutTag_Report.pdf
# - results/report/TES_TEAD1_CutTag_Report.pdf
#===============================================================================

#SBATCH --job-name=homer_to_pdf
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/homer_to_pdf.err"
#SBATCH --output="./logs/homer_to_pdf.out"

set -e

# Base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
# COMMENTED OUT: HOMER_DIR="${BASE_DIR}/results/08_homer_motifs"

# Activate conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if weasyprint environment exists, create if not
ENV_NAME="weasyprint_env"
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment with weasyprint..."
    conda create -n ${ENV_NAME} python=3.10 -y
    conda activate ${ENV_NAME}
    pip install weasyprint
else
    conda activate ${ENV_NAME}
fi

echo "=========================================="
echo "HTML to PDF Converter"
echo "Date: $(date)"
echo "=========================================="

# # COMMENTED OUT - Previous HOMER conditions
# # Conditions to process (can be overridden by command line argument)
# if [ -n "$1" ]; then
#     CONDITIONS=("$1")
# else
#     CONDITIONS=("TES" "TEAD1" "TESmut")
# fi

# Python script for conversion
python3 << 'PYTHON_SCRIPT'
import os
import sys
from pathlib import Path

# Try to import weasyprint
try:
    from weasyprint import HTML, CSS
    HAS_WEASYPRINT = True
except ImportError:
    HAS_WEASYPRINT = False
    print("Warning: weasyprint not available, trying alternative methods...")

def convert_with_weasyprint(html_path, pdf_path):
    """Convert HTML to PDF using weasyprint."""
    # Custom CSS for better PDF rendering
    custom_css = CSS(string='''
        @page {
            size: A4 landscape;
            margin: 1cm;
        }
        body {
            font-family: Arial, sans-serif;
            font-size: 10pt;
        }
        table {
            font-size: 8pt;
            border-collapse: collapse;
            width: 100%;
        }
        td, th {
            padding: 2px 4px;
            border: 1px solid #ccc;
        }
        h1 {
            font-size: 14pt;
            margin-bottom: 10px;
        }
        svg {
            max-width: 100%;
        }
    ''')

    try:
        html = HTML(filename=str(html_path))
        html.write_pdf(str(pdf_path), stylesheets=[custom_css])
        return True
    except Exception as e:
        print(f"  Error converting {html_path}: {e}")
        return False

def main():
    base_dir = Path('/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG')

    # HTML files to convert (relative to base_dir)
    html_files = [
        'results/multiqc/summary_report.html',
        'results/multiqc/TES_TEAD1_CutTag_Report.html',
        'results/report/TES_TEAD1_CutTag_Report.html',
    ]

    # # COMMENTED OUT - Previous HOMER file list
    # html_files = ['knownResults.html', 'homerResults.html']
    # conditions_str = os.environ.get('CONDITIONS', 'TES TEAD1 TESmut')
    # conditions = conditions_str.split()

    if not HAS_WEASYPRINT:
        print("ERROR: weasyprint is required but not installed.")
        print("Please install it with: pip install weasyprint")
        sys.exit(1)

    converted = 0
    failed = 0

    for html_file in html_files:
        html_path = base_dir / html_file
        pdf_path = html_path.with_suffix('.pdf')

        if not html_path.exists():
            print(f"\nSkipping {html_file}: file not found")
            continue

        print(f"\nConverting {html_file}...", end=' ')

        if convert_with_weasyprint(html_path, pdf_path):
            print(f"OK -> {pdf_path.name}")
            converted += 1
        else:
            print("FAILED")
            failed += 1

    # # COMMENTED OUT - Previous HOMER condition-based conversion
    # for condition in conditions:
    #     condition_dir = Path(base_dir) / f"{condition}_peaks"
    #
    #     if not condition_dir.exists():
    #         print(f"\nSkipping {condition}: directory not found")
    #         continue
    #
    #     print(f"\nProcessing {condition}...")
    #
    #     for html_file in html_files:
    #         html_path = condition_dir / html_file
    #         pdf_path = condition_dir / html_file.replace('.html', '.pdf')
    #
    #         if not html_path.exists():
    #             print(f"  Skipping {html_file}: file not found")
    #             continue
    #
    #         print(f"  Converting {html_file}...", end=' ')
    #
    #         if convert_with_weasyprint(html_path, pdf_path):
    #             print(f"OK -> {pdf_path.name}")
    #             converted += 1
    #         else:
    #             print("FAILED")
    #             failed += 1

    print(f"\n========================================")
    print(f"Conversion complete!")
    print(f"  Converted: {converted}")
    print(f"  Failed: {failed}")
    print(f"========================================")

if __name__ == '__main__':
    main()

PYTHON_SCRIPT

echo ""
echo "PDF conversion complete!"
echo "Date: $(date)"
