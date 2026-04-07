#!/usr/bin/env python3
"""
Export all supplementary tables as Excel sheets.
"""
import pandas as pd
import glob

def main():
    writer = pd.ExcelWriter("outputs/Supplementary_Tables.xlsx", engine="openpyxl")
    # List all CSV outputs
    csv_files = glob.glob("outputs/*.csv")
    for f in csv_files:
        sheet_name = f.split("/")[-1].replace(".csv", "")[:31]  # Excel sheet max 31 chars
        df = pd.read_csv(f)
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.close()
    print("Supplementary tables created.")

if __name__ == "__main__":
    main()