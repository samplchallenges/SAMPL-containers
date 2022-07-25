docker build -t robbason/score-coords:latest .
docker tag robbason/score-coords:latest ghcr.io/robbason/score-coords:latest
docker push ghcr.io/robbason/score-coords:latest
