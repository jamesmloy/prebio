FROM ubuntu:20.04 as build_env

RUN apt-get update && \
    apt-get install -y \
    build-essential ninja-build libblas-dev liblapack-dev \
    curl g++-10 wget

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py310_23.5.2-0-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV PATH /opt/conda/bin:$PATH

COPY requirements.yaml /prebio/requirements.yaml

WORKDIR /prebio/

RUN conda install -y -n base conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda env update --file requirements.yaml -n base

WORKDIR /deps

RUN git clone https://github.com/project-gemmi/gemmi.git && \
    mkdir gemmi/build && cd gemmi/build && \
    cmake -DPython_EXECUTABLE=$(which python) -DCMAKE_CXX_COMPILER=g++-10 -DUSE_PYTHON=ON .. && \
    make -j 5 &&  make install

FROM build_env as build_prebio

WORKDIR /build-prebio
COPY ./prebio ./prebio
COPY ./CMakeLists.txt ./
COPY ./pyproject.toml ./
COPY ./requirements.yaml ./
RUN pip3 wheel -w dist -vvv .

RUN conda install -y -c conda-forge conda-pack

RUN conda-pack -o /tmp/prod_env.tar.gz && \
    mkdir /venv && \
    cd /venv && \
    tar -xzf /tmp/prod_env.tar.gz && \
    rm /tmp/prod_env.tar.gz

RUN /venv/bin/conda-unpack

FROM ubuntu:20.04 as final

COPY --from=build_prebio /venv /venv
COPY --from=build_prebio /usr/local/include/. /usr/local/include
COPY --from=build_prebio /usr/local/lib/. /usr/local/lib
COPY --from=build_prebio /opt/conda/lib/python3.10/site-packages/gemmi.cpython-310-x86_64-linux-gnu.so /venv/lib/python3.10/site-packages

COPY --from=build_prebio /usr/local/bin/. /usr/local/bin
COPY --from=build_prebio /usr/local/share/. /usr/local/share
COPY ./generate_boxes_docker.py .
COPY --from=build_prebio /build-prebio/dist /prebio/dist

ENV PATH=/venv/bin:$PATH
ENV LD_LIBRARY_PATH /usr/local/lib:/venv/lib
RUN pip3 install /prebio/dist/$(ls /prebio/dist)

FROM final as dev

RUN apt-get update --fix-missing && apt-get install -y build-essential git \
    curl ninja-build g++-10 wget libblas-dev liblapack-dev && \
    curl -fsSL https://get.docker.com -o get-docker.sh && sh ./get-docker.sh

RUN pip uninstall -y prebio && \
    pip install pytest pytest-asyncio mypy scikit-build-core pybind11 cmake \
    pathspec packaging pyproject_metadata
ENV ZLIB_ROOT=/venv