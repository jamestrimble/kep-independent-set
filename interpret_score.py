import sys

x = long(sys.argv[1])
criteria = reversed(["effective",
                     "size     ",
                     "3way     ",
                     "backarc  ",
                     "weight   "])

for criterion in criteria:
    print criterion, (1023 - x % 1024) if criterion.startswith("3way") else x % 1024
    x /= 1024
