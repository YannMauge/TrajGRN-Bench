# Build all benchmark Docker images in one command (make -j4 all for parallel builds).
# All targets must be invoked from the repository root,
# e.g.  make -f containers/Makefile all
#
# Image definitions are sourced from methods_registry.yaml via
#   python utils/methods_registry.py container-images
# To add a new image, edit methods_registry.yaml — no Makefile changes needed.

REPO_ROOT  := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))..
DOCKERFILE := containers/docker/Dockerfile
DOCKER     := docker

# ── Image definitions (auto-generated from registry) ────────────────
define IMAGES_JSON
$(shell cd $(REPO_ROOT) && python utils/methods_registry.py container-images 2>/dev/null || echo '[]')
endef

# ── Extract target names ────────────────────────────────────────────
TARGETS := $(shell cd $(REPO_ROOT) && python utils/methods_registry.py container-images 2>/dev/null | python -c "import json,sys; [print(d['target']) for d in json.load(sys.stdin)]" 2>/dev/null)

# ── Phony targets ──────────────────────────────────────────────────
.PHONY: all list clean help $(TARGETS)

# ── all: build every image ─────────────────────────────────────────
all: $(TARGETS)

# ── Individual image targets (dynamically generated) ────────────────
define BUILD_RULE
$(1):
	@cd $(REPO_ROOT) && python -c "\
import json, subprocess, sys; \
images = json.loads(subprocess.check_output([sys.executable, 'utils/methods_registry.py', 'container-images'], text=True)); \
img = next((d for d in images if d['target'] == '$(1)'), None); \
if img: \
    subprocess.check_call(['$(DOCKER)', 'build', \
        '-f', '$(DOCKERFILE)', \
        '--build-arg', f\"ENV_FILE={img['env_file']}\", \
        '--build-arg', f\"ENV_NAME={img['env_name']}\", \
        '-t', f\"{img['image']}:latest\", '.']) \
"
endef

$(foreach t,$(TARGETS),$(eval $(call BUILD_RULE,$(t))))

# ── Utility targets ────────────────────────────────────────────────

list:  ## List all images that will be built
	@cd $(REPO_ROOT) && python -c "\
import json, subprocess, sys; \
images = json.loads(subprocess.check_output([sys.executable, 'utils/methods_registry.py', 'container-images'], text=True)); \
print('Images that will be built by all:'); \
[print(f\"  {d['image']}:latest  ({d['env_file']})\") for d in images]"

clean:  ## Remove all built benchmark images (if they exist)
	@cd $(REPO_ROOT) && python -c "\
import json, subprocess, sys; \
images = json.loads(subprocess.check_output([sys.executable, 'utils/methods_registry.py', 'container-images'], text=True)); \
for d in images: \
    img = d['image'] + ':latest'; \
    try: \
        subprocess.check_call(['$(DOCKER)', 'image', 'inspect', img], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL); \
        print(f'Removing {img}'); \
        subprocess.check_call(['$(DOCKER)', 'rmi', img]); \
    except: pass"

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-25s\033[0m %s\n", $$1, $$2}'
