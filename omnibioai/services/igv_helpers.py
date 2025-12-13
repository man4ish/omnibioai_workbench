import os
import json
import time
from selenium import webdriver
from selenium.webdriver.chrome.options import Options


def write_igv_html(
    output_html_path: str,
    genome: str = "hg38",
    tracks: list[dict] = None,
    locus: str = None
):
    """
    Generate an IGV.js HTML file for interactive genome visualization.

    Args:
        output_html_path (str): Path to save the HTML file.
        genome (str, optional): Genome assembly (default "hg38").
        tracks (list of dict, optional): List of IGV track dictionaries.
        locus (str, optional): Initial locus to display, e.g., "chr1:100000-200000".
    """
    tracks = tracks or []

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>IGV.js Genome Browser</title>
    <script src="https://cdn.jsdelivr.net/npm/igv/dist/igv.min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; }}
        #igvDiv {{ width: 100%; height: 600px; border: 1px solid black; }}
    </style>
</head>
<body>
    <div id="igvDiv"></div>
    <script>
        const options = {{
            genome: "{genome}",
            locus: "{locus or ''}",
            tracks: {json.dumps(tracks)}
        }};
        igv.createBrowser(document.getElementById("igvDiv"), options)
            .then(browser => {{
                console.log("IGV browser initialized:", browser);
            }});
    </script>
</body>
</html>
"""
    os.makedirs(os.path.dirname(output_html_path), exist_ok=True)
    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    print(f"[INFO] IGV HTML saved at {output_html_path}")


def capture_igv_snapshot(html_file: str, output_png: str, wait_time: int = 10):
    """
    Use headless Chrome to open an IGV HTML file and capture a PNG screenshot.

    Args:
        html_file (str): Path to IGV HTML file.
        output_png (str): Path to save screenshot PNG.
        wait_time (int, optional): Seconds to wait for IGV to load (default 10).
    """
    options = Options()
    options.headless = True
    options.add_argument("--window-size=1600,800")
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")

    driver = webdriver.Chrome(options=options)
    driver.get(f"file://{os.path.abspath(html_file)}")

    # Wait for IGV browser to load completely
    time.sleep(wait_time)

    # Capture screenshot
    driver.save_screenshot(output_png)
    driver.quit()
    print(f"[INFO] IGV snapshot saved at {output_png}")
