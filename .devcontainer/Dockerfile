# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.231.6/containers/cpp/.devcontainer/base.Dockerfile

# [Choice] Debian / Ubuntu version (use Debian 11, Ubuntu 18.04/21.04 on local arm64/Apple Silicon): debian-11, debian-10, ubuntu-21.04, ubuntu-20.04, ubuntu-18.04
ARG VARIANT="bullseye"
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

RUN apt update && \
	apt install -y \
		libeigen3-dev \
		doxygen && \
	rm -rf /var/lib/apt/lists/*

# Install Catch2
RUN git clone https://github.com/catchorg/Catch2.git && \
	cd Catch2 && \
	cmake -B build -S . -D BUILD_TESTING=OFF && \
	sudo cmake --build build/ --target install && \
	cd .. && \
	rm -r Catch2
