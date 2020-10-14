#!/usr/bin/env python 

import pandas as pd
import argparse as ap
from sqlalchemy import create_engine

# script to add tab delim file to sqlite database schema

parser = ap.ArgumentParser(description='add file to sqlite database')
parser.add_argument('a', help='path and file to add to database')
parser.add_argument('b', help='absolute path to database file')
parser.add_argument('c', help='Name of database table')
args = parser.parse_args()

filein = args.a


df = pd.read_csv(filein, header=0, sep="\t")

engine = create_engine('sqlite:///' + args.b)
with engine.connect() as conn, conn.begin():
    df.to_sql(args.c, conn, if_exists='append', index=False)
