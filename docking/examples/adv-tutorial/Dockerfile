FROM adv-base

WORKDIR /opt/app/

COPY autodock.py main.py ./

ENV PATH="/root/.local/bin:$PATH"

RUN ls -l /opt/app

ENTRYPOINT ["python3", "/opt/app/main.py"]
