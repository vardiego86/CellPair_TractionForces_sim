FROM ubuntu:bionic

LABEL maintainer="Naeem Muhammad <naeem.muhammad@kuleuven.be>"

USER root
WORKDIR /root

SHELL [ "/bin/bash", "-c" ]

ARG PYTHON_VERSION_TAG=3.6.7
ARG LINK_PYTHON_TO_PYTHON3=1

ARG CMAKE_VERSION=3.14.3

# Existing lsb_release causes issues with modern installations of Python3
# https://github.com/pypa/pip/issues/4924#issuecomment-435825490
# Set (temporarily) DEBIAN_FRONTEND to avoid interacting with tzdata
RUN apt-get -qq -y update && \
    apt-get -qq -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get -qq -y install \
        gcc=4:7.4.0-1ubuntu2.3 \
        g++=4:7.4.0-1ubuntu2.3 \
        zlibc \
        zlib1g-dev \
		texlive-latex-base \
		texlive-fonts-recommended \
		texlive-fonts-extra \
		texlive-latex-extra \
        libssl-dev \
        libbz2-dev \
        libsqlite3-dev \
        libncurses5-dev \
        libgdbm-dev \
        libgdbm-compat-dev \
        liblzma-dev \
        libreadline-dev \
        uuid-dev \
        libffi-dev \
        libtbb-dev \
        libcgal-dev \
        tk-dev \
        wget \
        curl \
        git \
        make \
        sudo \
        bash-completion \
        tree \
        vim \
        clang \
        libboost-dev=1.65.1.0ubuntu1 \
        libboost-all-dev=1.65.1.0ubuntu1 \
        #imagemagick=8:6.9.7.4+dfsg-16ubuntu6.7 \
        imagemagick \
        #shapely \
        texlive-extra-utils=2017.20180305-2 \
        python-h5py=2.7.1-2 \
        software-properties-common && \
    mv /usr/bin/lsb_release /usr/bin/lsb_release.bak && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/*

# Download cmake (ccmak) bash script and install
RUN cd /opt && wget https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
  && chmod +x cmake-${CMAKE_VERSION}-Linux-x86_64.sh \
  && bash ./cmake-${CMAKE_VERSION}-Linux-x86_64.sh < <(echo y) >/dev/null < <(echo y) \
  && ln -s /opt/cmake-${CMAKE_VERSION}-Linux-x86_64/bin/* /usr/bin \
&& rm -rf cmake-${CMAKE_VERSION}-Linux-x86_64.sh

COPY install_python.sh /tmp/install_python.sh
RUN bash /tmp/install_python.sh ${PYTHON_VERSION_TAG} ${LINK_PYTHON_TO_PYTHON3} && \
    rm -r /tmp/install_python.sh Python-${PYTHON_VERSION_TAG} && \
	ln -s /usr/bin/python3.6 /usr/bin/python-mpacts

RUN pip install Sphinx==2.2.0 vtk==8.1.2 numpy==1.17.2 networkx==2.3 scipy==1.3.1 matplotlib==2.2.4 h5py==2.10.0 pyevtk==1.1.1 seaborn==0.9.0 Shapely==1.6.4.post2 git+https://github.com/tonysyu/mpltools.git git+https://github.com/tisimst/pyDOE.git


# Enable tab completion by uncommenting it from /etc/bash.bashrc
# The relevant lines are those below the phrase "enable bash completion in interactive shells"
RUN export SED_RANGE="$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+1)),$(($(sed -n '\|enable bash completion in interactive shells|=' /etc/bash.bashrc)+7))" && \
    sed -i -e "${SED_RANGE}"' s/^#//' /etc/bash.bashrc && \
    unset SED_RANGE

# Create user "docker" with sudo powers
RUN useradd -m docker && \
    usermod -aG sudo docker && \
    echo '%sudo ALL=(ALL) NOPASSWD: ALL' >> /etc/sudoers && \
    cp /root/.bashrc /home/docker/ && \
    chown -R --from=root docker /home/docker

ENV \
    CC="/usr/bin/clang" \
    CXX="/usr/bin/clang++"

# Use C.UTF-8 locale to avoid issues with ASCII encoding
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8


WORKDIR /home/docker
ENV HOME /home/docker
ENV USER docker

RUN cd /tmp && wget -O mpacts-linux.deb --post-data "anonymPassword=mpacts" "https://drives.kuleuven.be/Handlers/AnonymousDownload.ashx?file=792b61c1&links=true" && \
	dpkg -i /tmp/mpacts-linux.deb && \
	rm /tmp/mpacts-linux.deb


USER docker
ENV PATH /home/docker/.local/bin:$PATH
# Avoid first use of sudo warning. c.f. https://askubuntu.com/a/22614/781671
RUN touch $HOME/.sudo_as_admin_successful


ENV PATH $HOME/bin/:$PATH
ENV PYTHONPATH $HOME/lib/:$PYTHONPATH

CMD [ "/bin/bash" ]
