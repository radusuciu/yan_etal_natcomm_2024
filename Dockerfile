ARG THERMO_RAWFILEPARSER_VERSION=1.4.3
ARG COMET_VERSION=2019015
ARG OPENMS_VERSION=3.1.0
ARG OPENMS_INSTALL_DIR="/opt/OpenMS"
ARG PERCOLATOR_VERSION=3-05
ARG R_VERSION=4.3.3
ARG DEBIAN_FRONTEND=noninteractive
ARG SAGE_VERSION=v0.14.3


### various software fetcher: comet etc. ### 
FROM debian:bullseye-slim AS downloader
ARG DEBIAN_FRONTEND
ARG THERMO_RAWFILEPARSER_VERSION
ARG COMET_VERSION
ARG SAGE_VERSION

RUN apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    # additional OS packages
    ca-certificates \
    wget \
    git \
    unzip \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace/thermorawfileparser_download
RUN wget https://github.com/compomics/ThermoRawFileParser/releases/download/v${THERMO_RAWFILEPARSER_VERSION}/ThermoRawFileParser${THERMO_RAWFILEPARSER_VERSION}.zip
RUN unzip *.zip && rm *.zip
RUN mv /workspace/thermorawfileparser_download/ /opt/ThermoRawFileParser/

WORKDIR /workspace/comet_download
RUN wget https://master.dl.sourceforge.net/project/comet-ms/comet_${COMET_VERSION}.zip
RUN unzip *.zip
RUN mv comet*debian* /usr/local/bin/comet && chmod +x /usr/local/bin/comet && rm -rf ../comet_download

WORKDIR /workspace/sage_download
RUN wget https://github.com/lazear/sage/releases/download/${SAGE_VERSION}/sage-${SAGE_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN tar -xzf *.tar.gz
RUN mv sage-${SAGE_VERSION}-x86_64-unknown-linux-gnu/sage /usr/local/bin/sage && chmod +x /usr/local/bin/sage && rm -rf ../sage_download


### percolator -- using a separate stage so we can use an arg to define the version ###
FROM ghcr.io/radusuciu/docker-percolator:${PERCOLATOR_VERSION} AS percolator-build


### openms -- using a separate stage so we can use an arg to define the version ###
FROM ghcr.io/radusuciu/docker-openms:${OPENMS_VERSION} AS openms


### worker ###
FROM python:3.11.9-slim-bullseye AS worker
ARG DEBIAN_FRONTEND
ARG R_VERSION
ARG OPENMS_INSTALL_DIR
ENV CIMAGE_PATH=/opt/cimage
ARG UID=1000
ARG GID=1000
ARG WORKER_USER=brain

ENV VIRTUAL_ENV=/venv
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8
ENV PATH="${OPENMS_INSTALL_DIR}/bin:${PATH}"

# create new user which will actually run the application
RUN addgroup --gid ${GID} ${WORKER_USER}
RUN adduser --disabled-password --gecos '' --uid ${UID} --gid ${GID} ${WORKER_USER}
RUN chown -R ${WORKER_USER} /home/${WORKER_USER}

RUN apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    gpg \
    gpg-agent \
    dirmngr \
  && sed -r -i 's/^deb(.*)$/deb\1 contrib/g' /etc/apt/sources.list \
  && echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections \
  && gpg --batch --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys B8F25A8A73EACF41 \
  && gpg --export --armor B8F25A8A73EACF41 | gpg --dearmor -o /usr/share/keyrings/marutter-archive-keyring.gpg \
  && echo 'APT::Sandbox::User "root";' > /etc/apt/apt.conf.d/99sandbox \
  && echo "deb [signed-by=/usr/share/keyrings/marutter-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/debian bullseye-cran40/" > /etc/apt/sources.list.d/cran.list \
  && echo "deb-src [signed-by=/usr/share/keyrings/marutter-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/debian bullseye-cran40/" >> /etc/apt/sources.list.d/cran.list \
  && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF \
  && echo "deb https://download.mono-project.com/repo/debian stable-buster main" > /etc/apt/sources.list.d/mono-official-stable.list \
  && apt-get update \
  && apt-get install -y --no-install-recommends --no-install-suggests \
    # OpenMS runtime deps
    libqt5opengl5 \
    libsvm3 \
    libzip4 \
    zlib1g \
    libbz2-1.0 \
    libgomp1 \
    libqt5svg5 \
    libxerces-c3.2 \
    coinor-libcoinmp1v5 \
    # cimage runtime deps
    r-base=${R_VERSION}* \
    libglpk40 \
    netcdf-bin \
    libxml2 \
    libnetcdf-dev \
    imagemagick \
    pdftk \
    # additional OS packages
    mono-complete \
    gnupg \
    gnupg2 \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    locales \
    locales-all \
    build-essential \
    git \
    unzip \
    rsync \
    ttf-mscorefonts-installer \
    fontconfig \
    libldap-2.4-2 \
    libsasl2-2 \
    libpoppler102 \
    libcairo2 \
    libpq5 \
    # for pyopenms
    libglib2.0-0 \
  && fc-cache -f \
  && rm -rf /var/lib/apt/lists/*

# cimage R dependencies
COPY --from=ghcr.io/radusuciu/docker-cimage-base /usr/local/lib/R/ /usr/local/lib/R/
COPY --from=ghcr.io/radusuciu/docker-cimage-base /usr/lib/R/ /usr/lib/R/
COPY --from=ghcr.io/radusuciu/docker-cimage-base /usr/share/R/ /usr/share/R/

# percolator
WORKDIR /tmp
COPY --from=percolator-build /workspace/release/*.deb .
RUN dpkg -i *.deb && rm -rf /tmp/*

# comet + ThermoRawFileParser + sage
COPY --from=downloader /usr/local/bin/* /usr/local/bin
COPY --from=downloader /opt/ThermoRawFileParser /opt/ThermoRawFileParser

# copy openms library and tools
COPY --from=openms ${OPENMS_INSTALL_DIR} ${OPENMS_INSTALL_DIR}
RUN chown -R ${UID}:${GID} ${OPENMS_INSTALL_DIR}/share

USER ${WORKER_USER}

# copy cimage files
COPY ./cimage ${CIMAGE_PATH}
LABEL org.opencontainers.image.source https://github.com/radusuciu/yan_etal_natcomm_2024
