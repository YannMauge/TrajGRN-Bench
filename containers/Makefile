# Thin wrapper around containers/build_all.sh (run from the repository root):
#   make -f containers/Makefile all          # all CPU images, native platform
#   make -f containers/Makefile all-cuda     # + CUDA images (linux/amd64)
#   make -f containers/Makefile flecs scnode # selected pixi envs
#   make -f containers/Makefile list
# Images: benchmark/<pixi_env>:latest, built from pixi.toml / pixi.lock.

BUILD := containers/build_all.sh
J ?= 1

.PHONY: all all-cuda list help

all:
	$(BUILD) -j $(J)

all-cuda:
	$(BUILD) -j $(J) --cuda

list:
	$(BUILD) --list

help:
	@sed -n '1,7p' containers/Makefile

%:
	$(BUILD) $@
