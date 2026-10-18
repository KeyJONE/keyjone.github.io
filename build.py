"""Build jemdoc sources and copy self-hosted assets into docs/ for GitHub Pages."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit
from generate_papers import generate as generate_papers

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "site"
TOOL = ROOT / "tools" / "jemdoc-cvx-1.0.0"
OUT = ROOT / "docs"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.paths = []
        self.invalid_targets = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ("href", "src") and value:
                parts = urlsplit(value)
                if not parts.scheme and not parts.netloc and parts.path:
                    self.paths.append(unquote(parts.path))
            if key == "target" and value != "_blank":
                self.invalid_targets.append(value)


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
    generate_papers()
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
    paper_sources = sorted((SITE / "papers").glob("*.jemdoc"))
    if paper_sources:
        (OUT / "papers").mkdir(exist_ok=True)
        subprocess.run(
            [sys.executable, "-X", "utf8", str(TOOL / "jemdoc"), "-c", "../papers.conf",
             "-o", str(OUT / "papers") + os.sep, *[p.name for p in paper_sources]],
            cwd=SITE / "papers", env=env, check=True,
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
            shutil.copytree(SITE / directory, OUT / directory, dirs_exist_ok=True,
                            ignore=shutil.ignore_patterns("*.jemdoc", "*.conf", "katexcache", "__pycache__"))
    (OUT / ".nojekyll").touch()
    all_sources = [SITE / source for source in sources] + paper_sources
    for source in all_sources:
        page = OUT / source.relative_to(SITE).with_suffix(".html")
        markup = page.read_text(encoding="utf-8")
        if 'class="katex-error"' in markup:
            raise RuntimeError(f"Invalid LaTeX in {source}; please fix it before publishing.")
        original = source.read_text(encoding="utf-8")
        if ("$" in original or "\\(" in original) and 'class="katex"' not in markup:
            raise RuntimeError(f"No rendered math found in {source}; check KaTeX.")
        links = Links()
        links.feed(markup)
        if links.invalid_targets:
            raise RuntimeError(f"Invalid link target in {source}: {links.invalid_targets}")
        for relative in links.paths:
            if not (page.parent / relative).exists():
                raise RuntimeError(f"Broken local link in {source}: {relative}")
    print(f"Built and checked {len(all_sources)} pages. Open: {OUT / 'index.html'}")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"Build failed: {exc}", file=sys.stderr)
        sys.exit(1)
