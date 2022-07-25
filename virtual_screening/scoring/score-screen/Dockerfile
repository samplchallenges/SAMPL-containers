FROM ghcr.io/megosato/basepyrd:latest

WORKDIR /opt/app

COPY score.py ./

ENV PATH="/root/.local/bin:$PATH"

ENTRYPOINT ["python3", "/opt/app/score.py"]
