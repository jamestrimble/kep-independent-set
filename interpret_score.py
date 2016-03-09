import sys

x = long(sys.argv[1])
criteria = reversed(["effective   ",
                     "size        ",
                     "inverse3way ",
                     "backarc     ",
                     "weight      "])

for criterion in criteria:
    print criterion, x % 1024
    x /= 1024
