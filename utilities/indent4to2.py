#!/usr/bin/env python3

import sys
import re
from pathlib import Path

if len(sys.argv) != 2:
  print(f"Usage: {sys.argv[0]} FILE")
  sys.exit(1)

path = Path(sys.argv[1])

if not path.is_file():
  print(f"Error: '{path}' is not a file")
  sys.exit(1)

text = path.read_text()

new_lines = []

for line in text.splitlines(keepends=True):
  match = re.match(r"^( +)", line)

  if match:
    n = len(match.group(1))

    # Convert indentation levels of 4 spaces to 2 spaces
    if n % 4 == 0:
      line = " " * (n // 2) + line[n:]

  new_lines.append(line)

path.write_text("".join(new_lines))
