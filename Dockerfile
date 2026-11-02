# syntax=docker/dockerfile:1

# ---- toolchain: g++, cmake and Eigen, shared by the C++ and Python builds ----
FROM ubuntu:24.04 AS toolchain

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src

# ---- build: compile the library, app and tests ----
FROM toolchain AS build

COPY CMakeLists.txt ./
COPY include/ include/
COPY src/ src/
COPY app/ app/
COPY tests/ tests/

# the tests rely on assert(), so keep it enabled by dropping Release's default -DNDEBUG
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE="-O2" \
    && cmake --build build --parallel

# ---- test: `docker build --target test .` fails if any unit test fails ----
FROM build AS test
RUN ctest --test-dir build --output-on-failure

# ---- python: `docker run --rm -it <image>` opens Python with spinchain installed ----
FROM toolchain AS python

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.13.6 /uv /usr/local/bin/uv

COPY CMakeLists.txt pyproject.toml uv.lock README.md ./
COPY include/ include/
COPY src/ src/
COPY python/ python/

RUN uv sync --locked --no-editable --python /usr/bin/python3

CMD [".venv/bin/python"]

# ---- test-python: `docker build --target test-python .` fails if any test fails ----
FROM python AS test-python
RUN .venv/bin/python -m pytest

# ---- runtime: just the app executable ----
FROM ubuntu:24.04 AS runtime

# the app links against libstdc++ dynamically
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/build/app/main /usr/local/bin/main

ENTRYPOINT ["main"]
