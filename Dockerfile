FROM rocker/r-base
LABEL maintainer="Vy Nguyen <vn@genomicexpression.com>"
ENV DEBIAN_FRONTEND noninteractive
ARG R_VERSION=4.1.3

COPY . /app
RUN mkdir -p /root/.aws
COPY /credentials /root/.aws
COPY /config /root/.aws

RUN apt-get update -y && apt-get install -y \
    curl \
    python3-dev \
    python3-boto3 \
    && mkdir -p /software

RUN apt-get update && \
 apt-get install -y \
    nodejs npm

WORKDIR /app/lib/src
RUN npm install csv-parse

COPY /lib /software

## Clean up
RUN cd / \
    && rm -rf /tmp/* \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean \
    && rm -rf /app/config \
    && rm -rf /app/credentials

# update path
ENV PATH="/app/:${PATH}"

CMD cd ~/../app && python3 handler.py

