from pathlib import Path
import sys
import pypdfium2 as pdfium
pdf = pdfium.PdfDocument(str(Path(sys.argv[1])))
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
for index in range(len(pdf)):
    page = pdf[index]
    image = page.render(scale=2.0).to_pil()
    image.save(out / f"page-{index + 1:03d}.png")
    page.close()
print(len(pdf))
