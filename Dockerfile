# syntax=docker/dockerfile:1

# ---- build: compile the library, app and tests ----
FROM ubuntu:24.04 AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
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

# ---- runtime: just the app executable ----
FROM ubuntu:24.04 AS runtime

# the app links against libstdc++ dynamically
RUN apt-get update && apt-get install -y --no-install-recommends \
        libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /src/build/app/main /usr/local/bin/main

ENTRYPOINT ["main"]
