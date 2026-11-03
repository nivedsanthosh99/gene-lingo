# Hugging Face Space (SDK: Docker). Free CPU Basic: 2 vCPU / 16 GB. See deploy/huggingface/README.md.
FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends diamond-aligner \
    && rm -rf /var/lib/apt/lists/*

RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user PATH=/home/user/.local/bin:$PATH
WORKDIR /home/user/app

COPY --chown=user . .
RUN pip install --no-cache-dir ".[ui]"

EXPOSE 7860
CMD ["streamlit", "run", "app.py", "--server.port=7860", "--server.address=0.0.0.0", "--server.headless=true"]
