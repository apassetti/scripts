#!/usr/bin/env python

import pdfplumber
import tkinter as tk
from tkinter import filedialog, messagebox
from transformers import pipeline
from fpdf import FPDF
import os

# Hide the root window and ask for a PDF file
root = tk.Tk()
root.withdraw()
pdf_file = filedialog.askopenfilename(filetypes=[("PDF files", "*.pdf")])
root.destroy()  # Close root window

if not pdf_file:
    print("No file selected.")
else:
    try:
        # Open the PDF and extract text from all pages
        with pdfplumber.open(pdf_file) as pdf:
            text = ""
            for page in pdf.pages:
                page_text = page.extract_text()
                if page_text:
                    text += page_text + "\n"
        
        # Check if the folder exists, then create if it doesn't
        output_dir = "/Users/andreapassetti/tmp/Supramolecular"
        os.makedirs(output_dir, exist_ok=True)
        
        # Save extracted text to a new PDF with Unicode support in ~/tmp folder
        output_pdf_path = os.path.join(output_dir, "summary_output.pdf")
        pdf_writer = FPDF()
        pdf_writer.add_page()
        pdf_writer.set_auto_page_break(auto=True, margin=15)
        pdf_writer.set_font("Helvetica", size=12)  # Switch to a core font

        # Define max width for multi_cell to handle wide text
        max_width = 180  # Adjust this width to fit page margins

        # Write text to PDF, handling special characters
        for line in text.split('\n'):
            pdf_writer.multi_cell(max_width, 10, line.encode('latin-1', 'replace').decode('latin-1'))

        pdf_writer.output(output_pdf_path)
        print(f"Extracted text saved to {output_pdf_path}")

        # Validate text length before summarizing
        if len(text.strip()) < 50:
            raise ValueError("Insufficient text found in PDF for summarization.")

        # Initialize the summarization pipeline, using GPU if available
        summarizer = pipeline("summarization", model="facebook/bart-large-cnn", device=0)

        # Summarize text (adjust max/min length as needed)
        summary = summarizer(text, max_length=2000, min_length=500, do_sample=False)
        print("\nSummary:")
        print(summary[0]['summary_text'])

    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
