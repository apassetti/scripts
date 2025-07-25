# Voice Notes - PDF to Audio Converter

This application converts PDF documents to audio files using:
- **PDF text extraction** with pdfplumber
- **Text processing** with Ollama (llama3.1 model) for grammar improvement
- **Text-to-speech** with Microsoft Edge TTS
- **Audio processing** with ffmpeg for long audio files

## Features

- Extract text from PDF files
- Improve text grammar and readability using AI
- Convert text to high-quality speech
- **Handles long documents** by splitting into chunks and combining audio
- Real-time progress tracking
- Customizable voice settings

## Installation

1. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Install system dependencies:
   ```bash
   sudo apt install ffmpeg  # For Linux
   # or brew install ffmpeg  # For macOS
   ```

3. Make sure Ollama is installed and the llama3.1 model is available

## Usage

Run the application:
```bash
python voicenotes.py
```

1. Click "Select PDF and Start"
2. Choose your PDF file
3. Wait for text extraction and processing
4. Choose where to save the MP3 file
5. The app will generate audio, handling long texts automatically

## Configuration

You can modify these settings in the code:
- **Voice**: Change `voice="en-US-AriaNeural"` to use different voices
- **Speed**: Adjust `rate="+20%"` to change speech speed
- **Chunk size**: Modify `max_chars=8000` in `split_text_for_audio()` for different audio segment lengths