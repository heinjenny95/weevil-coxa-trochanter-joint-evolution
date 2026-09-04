from __future__ import annotations

import argparse
from io import BytesIO
import os
from pathlib import Path
import re

from PIL import Image, ImageDraw, ImageFont
from pypdf import PdfReader, PdfWriter
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


LAYOUTS = {
    "main_figure_2": {
        "reference_size": (4250.0, 2037.0),
        "old_boxes": ((264, 60, 318, 120), (2347, 50, 2402, 120)),
        "new_positions": ((103, 66), (2184, 56)),
        "labels": ("a", "b"),
        "font_pixels": 44,
        "font_points": 6.5,
    },
    "supplementary_figure_2": {
        "reference_size": (3952.0, 4324.0),
        "old_boxes": (
            (102, 52, 154, 113),
            (2198, 43, 2253, 114),
            (102, 1448, 154, 1512),
            (2197, 1439, 2252, 1514),
            (102, 2844, 154, 2910),
            (2195, 2833, 2245, 2915),
        ),
        "new_positions": (
            (86, 61),
            (2181, 52),
            (58, 1457),
            (2153, 1448),
            (86, 2853),
            (2181, 2842),
        ),
        "labels": tuple("abcdef"),
        "font_pixels": 44,
        "font_points": 6.5,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Move the two legacy Figure 2 panel tags to the outer y-axis-title edge."
    )
    parser.add_argument("--input-png", type=Path, required=True)
    parser.add_argument("--input-pdf", type=Path, required=True)
    parser.add_argument("--output-png", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--input-svg", type=Path)
    parser.add_argument("--output-svg", type=Path)
    parser.add_argument(
        "--font-path",
        type=Path,
        help=(
            "Bold TrueType font used for PNG and PDF panel labels. If omitted, "
            "WEV_FIGURE_BOLD_FONT and common Arial/DejaVu Sans locations are checked."
        ),
    )
    parser.add_argument("--layout", choices=sorted(LAYOUTS), required=True)
    return parser.parse_args()


def resolve_font_path(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    env_font = os.environ.get("WEV_FIGURE_BOLD_FONT")
    if env_font:
        candidates.append(Path(env_font))
    windows_dir = os.environ.get("WINDIR")
    if windows_dir:
        candidates.append(Path(windows_dir) / "Fonts" / "arialbd.ttf")
    candidates.extend(
        [
            Path("/Library/Fonts/Arial Bold.ttf"),
            Path("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
        ]
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    checked = "\n  ".join(str(path) for path in candidates)
    raise FileNotFoundError(
        "No bold TrueType font was found. Supply --font-path or set "
        f"WEV_FIGURE_BOLD_FONT. Checked:\n  {checked}"
    )


def reposition_svg(source: Path, destination: Path, layout: dict[str, object]) -> None:
    svg = source.read_text(encoding="utf-8")
    view_box_match = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', svg)
    if view_box_match is None:
        raise ValueError("SVG does not contain the expected zero-origin viewBox")
    view_width = float(view_box_match.group(1))
    reference_width, _ = layout["reference_size"]
    for label, (x, _) in zip(layout["labels"], layout["new_positions"], strict=True):
        target_x = x / reference_width * view_width
        pattern = rf'(<text\b[^>]*\bx=")[0-9.]+("[^>]*>{re.escape(label)}</text>)'
        svg, count = re.subn(pattern, rf'\g<1>{target_x:.6f}\g<2>', svg, count=1)
        if count != 1:
            raise ValueError(f"Expected exactly one SVG panel tag {label!r}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(svg, encoding="utf-8")


def reposition_png(
    source: Path, destination: Path, layout: dict[str, object], font_path: Path
) -> None:
    image = Image.open(source).convert("RGB")
    draw = ImageDraw.Draw(image)
    reference_width, reference_height = layout["reference_size"]
    sx = image.width / reference_width
    sy = image.height / reference_height
    for x0, y0, x1, y1 in layout["old_boxes"]:
        draw.rectangle((x0 * sx, y0 * sy, x1 * sx, y1 * sy), fill="white")
    font = ImageFont.truetype(str(font_path), max(9, round(layout["font_pixels"] * sy)))
    for label, (x, y) in zip(layout["labels"], layout["new_positions"], strict=True):
        draw.text((x * sx, y * sy), label, fill="#263746", font=font)
    destination.parent.mkdir(parents=True, exist_ok=True)
    image.save(destination, dpi=image.info.get("dpi", (300, 300)))


def reposition_pdf(
    source: Path, destination: Path, layout: dict[str, object], font_path: Path
) -> None:
    reader = PdfReader(str(source))
    if len(reader.pages) != 1:
        raise ValueError("Expected a one-page Figure 2 PDF")
    page = reader.pages[0]
    width = float(page.mediabox.width)
    height = float(page.mediabox.height)
    reference_width, reference_height = layout["reference_size"]
    sx = width / reference_width
    sy = height / reference_height

    packet = BytesIO()
    overlay = canvas.Canvas(packet, pagesize=(width, height))
    pdfmetrics.registerFont(TTFont("Figure-Bold", str(font_path)))
    overlay.setFillColorRGB(1, 1, 1)
    for x0, y0, x1, y1 in layout["old_boxes"]:
        overlay.rect(x0 * sx, height - y1 * sy, (x1 - x0) * sx, (y1 - y0) * sy, fill=1, stroke=0)
    overlay.setFillColorRGB(38 / 255, 55 / 255, 70 / 255)
    overlay.setFont("Figure-Bold", layout["font_points"])
    baseline_offset = 35
    for label, (x, y) in zip(layout["labels"], layout["new_positions"], strict=True):
        overlay.drawString(x * sx, height - (y + baseline_offset) * sy, label)
    overlay.save()
    packet.seek(0)

    page.merge_page(PdfReader(packet).pages[0])
    writer = PdfWriter()
    writer.add_page(page)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("wb") as handle:
        writer.write(handle)


def main() -> None:
    args = parse_args()
    layout = LAYOUTS[args.layout]
    font_path = resolve_font_path(args.font_path)
    reposition_png(args.input_png, args.output_png, layout, font_path)
    reposition_pdf(args.input_pdf, args.output_pdf, layout, font_path)
    if (args.input_svg is None) != (args.output_svg is None):
        raise ValueError("--input-svg and --output-svg must be supplied together")
    if args.input_svg is not None:
        reposition_svg(args.input_svg, args.output_svg, layout)


if __name__ == "__main__":
    main()
