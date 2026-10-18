"""Build jemdoc sources and copy self-hosted assets into docs/ for GitHub Pages."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "site"
TOOL = ROOT / "tools" / "jemdoc-cvx-1.0.0"
OUT = ROOT / "docs"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ("href", "src") and value:
                parts = urlsplit(value)
                if not parts.scheme and not parts.netloc and parts.path:
                    self.paths.append(unquote(parts.path))


def main():
    if sys.version_info < (3, 11):
        raise RuntimeError("Please use Python 3.11 or newer.")
    runtime = ROOT / "runtime.local.json"
    settings = json.loads(runtime.read_text(encoding="utf-8")) if runtime.exists() else {}
    node = shutil.which("node") or settings.get("node")
    if not node or not Path(node).is_file():
        raise RuntimeError("Node.js is missing. Install Node.js, then build again.")
    katex = TOOL / "node_modules" / "katex"
    if not (katex / "cli.js").is_file():
        raise RuntimeError("KaTeX is missing. Run npm install inside tools/jemdoc-cvx-1.0.0 first.")
    sources = sorted(p.name for p in SITE.glob("*.jemdoc"))
    if "index.jemdoc" not in sources:
        raise RuntimeError("site/index.jemdoc is required.")
    OUT.mkdir(exist_ok=True)
    env = os.environ.copy()
    env["PATH"] = str(Path(node).parent) + os.pathsep + env.get("PATH", "")
    env["PYTHONUTF8"] = "1"
    subprocess.run(
        [sys.executable, "-X", "utf8", str(TOOL / "jemdoc"), "-c", "site.conf",
         "-o", str(OUT) + os.sep, *sources],
        cwd=SITE, env=env, check=True,
    )
    shutil.copy2(TOOL / "css" / "jemdoc-cvx.css", OUT)
    shutil.copytree(TOOL / "css" / "fonts", OUT / "fonts", dirs_exist_ok=True)
    (OUT / "katex").mkdir(exist_ok=True)
    shutil.copy2(katex / "dist" / "katex.min.css", OUT / "katex")
    shutil.copy2(katex / "LICENSE", OUT / "katex" / "LICENSE.txt")
    shutil.copytree(katex / "dist" / "fonts", OUT / "katex" / "fonts", dirs_exist_ok=True)
    for css in SITE.glob("*.css"):
        shutil.copy2(css, OUT)
    for directory in ("images", "papers", "files"):
        if (SITE / directory).is_dir():
            shutil.copytree(SITE / directory, OUT / directory, dirs_exist_ok=True)
    (OUT / ".nojekyll").touch()
    for source in sources:
        page = OUT / Path(source).with_suffix(".html")
        markup = page.read_text(encoding="utf-8")
        if 'class="katex-error"' in markup:
            raise RuntimeError(f"Invalid LaTeX in {source}; please fix it before publishing.")
        original = (SITE / source).read_text(encoding="utf-8")
        if ("$" in original or "\\(" in original) and 'class="katex"' not in markup:
            raise RuntimeError(f"No rendered math found in {source}; check KaTeX.")
        links = Links()
        links.feed(markup)
        for relative in links.paths:
            if not (OUT / relative).exists():
                raise RuntimeError(f"Broken local link in {source}: {relative}")
    print(f"Built and checked {len(sources)} pages. Open: {OUT / 'index.html'}")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"Build failed: {exc}", file=sys.stderr)
        sys.exit(1)
