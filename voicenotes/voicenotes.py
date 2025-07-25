#!/usr/bin/env python3

import pdfplumber
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from ollama import chat
import threading
import asyncio
import edge_tts

# -----------------------
# Helper Functions
# -----------------------

def extract_pdf_text(pdf_path, progress_label, progress_bar):
    text = ""
    try:
        with pdfplumber.open(pdf_path) as pdf:
            total_pages = len(pdf.pages)
            progress_bar.config(maximum=total_pages)
            for i, page in enumerate(pdf.pages, start=1):
                text += (page.extract_text() or "") + "\n"
                _update_gui(progress_label, f"Extracting page {i}/{total_pages}...")
                _update_gui(progress_bar, i, is_progress=True)
        return text
    except Exception as e:
        messagebox.showerror("Error", f"Failed to read PDF file: {e}")
        return None

def _update_gui(widget, text_or_value, is_progress=False):
    if is_progress:
        widget.after(0, lambda: widget.config(value=text_or_value))
    else:
        widget.after(0, lambda: widget.config(text=text_or_value))

def _append_output(output_area, content):
    output_area.after(0, lambda: (output_area.insert(tk.END, content),
                                  output_area.see(tk.END)))

async def synthesize_speech(text, output_file, rate="+20%", voice="en-US-AriaNeural"):
    """Synthesize speech in chunks to avoid 10-minute limit and combine them using ffmpeg."""
    import tempfile
    import os
    import subprocess
    
    # Split text into smaller chunks for audio generation (roughly 5-7 minutes each)
    audio_chunks = split_text_for_audio(text, max_chars=8000)
    
    if len(audio_chunks) == 1:
        # If only one chunk, use the original method
        communicate = edge_tts.Communicate(text=text, voice=voice, rate=rate)
        await communicate.save(output_file)
        return
    
    # Generate audio for each chunk
    temp_files = []
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Generate audio files for each chunk
        for i, chunk in enumerate(audio_chunks):
            temp_file = os.path.join(temp_dir, f"chunk_{i:03d}.mp3")
            communicate = edge_tts.Communicate(text=chunk, voice=voice, rate=rate)
            await communicate.save(temp_file)
            temp_files.append(temp_file)
        
        # Create a file list for ffmpeg
        filelist_path = os.path.join(temp_dir, "filelist.txt")
        with open(filelist_path, 'w') as f:
            for temp_file in temp_files:
                f.write(f"file '{temp_file}'\n")
        
        # Use ffmpeg to concatenate the audio files
        cmd = [
            'ffmpeg', '-f', 'concat', '-safe', '0', 
            '-i', filelist_path, '-c', 'copy', 
            '-y', output_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise Exception(f"FFmpeg failed: {result.stderr}")
        
    finally:
        # Clean up temporary files
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        if os.path.exists(filelist_path):
            os.remove(filelist_path)
        if os.path.exists(temp_dir):
            os.rmdir(temp_dir)

def split_text_for_audio(text, max_chars=8000):
    """Split text into chunks suitable for audio generation."""
    if len(text) <= max_chars:
        return [text]
    
    chunks = []
    sentences = text.split('. ')
    current_chunk = ""
    
    for sentence in sentences:
        # Add the sentence to current chunk if it fits
        if len(current_chunk) + len(sentence) + 2 <= max_chars:  # +2 for '. '
            if current_chunk:
                current_chunk += ". " + sentence
            else:
                current_chunk = sentence
        else:
            # If current chunk is not empty, save it and start a new one
            if current_chunk:
                chunks.append(current_chunk.strip())
                current_chunk = sentence
            else:
                # If single sentence is too long, split it by words
                words = sentence.split()
                for word in words:
                    if len(current_chunk) + len(word) + 1 <= max_chars:
                        if current_chunk:
                            current_chunk += " " + word
                        else:
                            current_chunk = word
                    else:
                        if current_chunk:
                            chunks.append(current_chunk.strip())
                            current_chunk = word
                        else:
                            # Single word is too long, just add it
                            chunks.append(word)
                            current_chunk = ""
    
    # Add the last chunk if it's not empty
    if current_chunk:
        chunks.append(current_chunk.strip())
    
    return chunks

def split_text_into_chunks(text, max_chunk_size=2000):
    """Split text into chunks while trying to preserve sentence boundaries."""
    chunks = []
    sentences = text.split('. ')
    current_chunk = ""
    
    for sentence in sentences:
        # Add the sentence to current chunk if it fits
        if len(current_chunk) + len(sentence) + 2 <= max_chunk_size:  # +2 for '. '
            if current_chunk:
                current_chunk += ". " + sentence
            else:
                current_chunk = sentence
        else:
            # If current chunk is not empty, save it and start a new one
            if current_chunk:
                chunks.append(current_chunk.strip())
                current_chunk = sentence
            else:
                # If single sentence is too long, split it by words
                words = sentence.split()
                for word in words:
                    if len(current_chunk) + len(word) + 1 <= max_chunk_size:
                        if current_chunk:
                            current_chunk += " " + word
                        else:
                            current_chunk = word
                    else:
                        if current_chunk:
                            chunks.append(current_chunk.strip())
                            current_chunk = word
                        else:
                            # Single word is too long, just add it
                            chunks.append(word)
                            current_chunk = ""
    
    # Add the last chunk if it's not empty
    if current_chunk:
        chunks.append(current_chunk.strip())
    
    return chunks

def process_with_ollama(text, output_area, progress_label, progress_bar):
    try:
        # Split text into chunks
        chunks = split_text_into_chunks(text, max_chunk_size=2000)
        total_chunks = len(chunks)
        
        _update_gui(progress_label, f"Processing text in {total_chunks} chunks with Ollama...")
        progress_bar.config(mode='determinate', maximum=total_chunks)
        
        processed_chunks = []
        
        # Process each chunk
        for i, chunk in enumerate(chunks, 1):
            _update_gui(progress_label, f"Processing chunk {i}/{total_chunks}...")
            _update_gui(progress_bar, i, is_progress=True)
            
            stream = chat(
                model='llama3.1',
                messages=[{
                    'role': 'user',
                    'content': 'Improve the grammar of the text, making it more readable. Add punctuation. Output directly the improved text. Do not include any explanations or additional text. Keep the text in english:\n\n' + chunk
                }],
                stream=True,
            )

            processed_chunk = ""
            for chunk_response in stream:
                content = chunk_response['message']['content']
                processed_chunk += content
                _append_output(output_area, content)
            
            processed_chunks.append(processed_chunk.strip())

        # Combine all processed chunks
        final_processed_text = " ".join(processed_chunks)
        
        _update_gui(progress_bar, total_chunks, is_progress=True)

        if final_processed_text.strip():
            output_file = filedialog.asksaveasfilename(
                defaultextension=".mp3",
                filetypes=[("MP3 files", "*.mp3")]
            )

            if output_file:
                _update_gui(progress_label, "Generating audio from combined text...")
                
                asyncio.run(synthesize_speech(
                    text=final_processed_text,
                    output_file=output_file,
                    rate="+20%",  # adjust speed here
                    voice="en-US-AriaNeural"  # change voice here
                ))

                _update_gui(progress_label, "Audio generation completed.")
                messagebox.showinfo("Success", f"Audio saved to {output_file}")
            else:
                _update_gui(progress_label, "Text processing completed (no audio saved).")

        return final_processed_text

    except Exception as e:
        progress_bar.stop()
        messagebox.showerror("Error", f"Processing failed: {e}")
        return None

def start_processing():
    pdf_file = filedialog.askopenfilename(filetypes=[("PDF files", "*.pdf")])
    if not pdf_file:
        messagebox.showerror("Error", "No PDF file selected.")
        return

    progress_label.config(text="Starting PDF extraction...")
    progress_bar['value'] = 0
    output_area.delete("1.0", tk.END)

    def task():
        text = extract_pdf_text(pdf_file, progress_label, progress_bar)
        if text:
            process_with_ollama(text, output_area, progress_label, progress_bar)

    threading.Thread(target=task, daemon=True).start()

# -----------------------
# GUI Setup
# -----------------------

root = tk.Tk()
root.title("PDF to Text with Audio (Ollama + Edge TTS)")

frame = tk.Frame(root, padx=10, pady=10)
frame.pack(fill="both", expand=True)

progress_label = tk.Label(frame, text="Select a PDF to start.", anchor="w")
progress_label.pack(fill="x", pady=(0, 5))

progress_bar = ttk.Progressbar(frame, orient="horizontal", length=400, mode="determinate")
progress_bar.pack(fill="x", pady=(0, 10))

output_area = scrolledtext.ScrolledText(frame, wrap=tk.WORD, height=20, width=80)
output_area.pack(fill="both", expand=True, pady=(0, 5))

button_frame = tk.Frame(frame)
button_frame.pack(fill="x")

start_button = tk.Button(button_frame, text="Select PDF and Start", command=start_processing)
start_button.pack(side="left", padx=5)

root.mainloop()
