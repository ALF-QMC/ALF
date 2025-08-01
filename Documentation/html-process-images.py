#!/usr/bin/env python3

from bs4 import BeautifulSoup
import subprocess
from glob import glob

with open("doc.html", "r", encoding="utf-8") as f:
    html = f.read()

soup = BeautifulSoup(html, "html.parser")

imgs = []

for tag in soup.find_all("embed", src=True):
    if tag["src"].lower().endswith(".pdf"):
        imgs.append(tag["src"])

for img in imgs:
    files = glob(f'**/{img}', recursive=True)
    if len(files) != 1:
        print(f'Could not find {img}: {files}')
        continue
    subprocess.run(["pdftocairo", "-png", "-singlefile", f'../{files[0]}'], cwd='html')
    html = html.replace(img, img.split('/')[-1].replace('.pdf', '.png'))

subprocess.run(["pdftocairo", "-png", "-singlefile", '../Figures/U_Profile.pdf'], cwd='html')
html = html.replace('U_profile.pdf', 'U_Profile.png')


with open("html/doc.html", "w", encoding="utf-8") as f:
    f.write(html)
