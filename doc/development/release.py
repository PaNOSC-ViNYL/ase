"""Translate user names to real names"""
import csv
import sys
names = {}
for name, uname, email in csv.reader(open('committers.csv')):
    names[uname] = name
allnames = []
for uname in sys.stdin.read().splitlines():
    allnames.append(names.get(uname, uname))
print(', '.join(sorted(allnames)))
