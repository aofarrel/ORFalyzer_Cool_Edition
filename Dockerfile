FROM python:3.9.10-slim

COPY orfalyzer/orfalyzer.py ./orfalyzer/
COPY orfalyzer/sequenceanalysis.py ./orfalyzer/
COPY mini_e_coli.fna .

# needed for the utf-8 output 
ENV LANG C.UTF-8

CMD [ "ls" ]