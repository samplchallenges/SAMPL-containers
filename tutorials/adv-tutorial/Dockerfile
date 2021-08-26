#FROM adv-rdkit-base:latest
FROM adv-tutorial-base

WORKDIR /opt/app/

ENV PATH="/root/.local/bin:$PATH"

COPY setup.py autodock.py main.py ./

RUN pip install .

RUN ls -l /opt/app

ENTRYPOINT ["run-autodock"]
