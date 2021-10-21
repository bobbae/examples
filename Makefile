.POSIX:
DESTDIR=public

#alias si := system-info

#https://github.com/casey/just#installation

#system-info:
#	@echo from ~/justfile
#	@echo "This is an {{arch()}} machine running {{os()}}".
#	@echo PATH is {{env_var("PATH")}}
#	@echo PATH is $PATH
#	@echo current invocation directory is {{invocation_directory()}}

.PHONY: all
all:
	@echo make all is all we do
.PHONY: get_repository

get_repository:
	@echo "🛎 Getting Pages repository"
	#git clone https://github.com/bobbae/examples.git $(DESTDIR)

.PHONY: clean
clean:
	@echo "🧹 Cleaning old build"
	#cd $(DESTDIR) && rm -rf *

.PHONY: get
get:
	@echo "❓ Checking for make"
	#@if ! [ -x "$$(command -v make)" ]; then\
	#	echo "🤵 Getting make";\
	#	sudo apt install make;\
	#fi

.PHONY: build
build:
	@echo "🍳 building"
	# do the building
	# @echo "🧂 Optimizing images"
	# do the optimizing


.PHONY: test
test:
	@echo "🍜 Testing"
	# do the testing
	#docker run -v $(GITHUB_WORKSPACE)/$(DESTDIR)/:/mnt bobbae/examples mnt --disable-external

.PHONY: deploy
deploy:
	@echo "🎁 Preparing commit"
	#@cd $(DESTDIR) \
	#	&& git config user.email "bobbae+githubaction@gmail.com" \
	#	&& git config user.name "My Examples via GitHub Actions" \
	#	&& git add . \
	#	&& git status \
	#	&& git commit -m "🤖 CD bot is helping" \
	#	&& git push -f -q https://$(TOKEN)@github.com/bobbae/examples master
	@echo "🚀 deployed!"

.PHONY: run
run:
	@echo "Run"
	# run
	@echo "Done running"

.PHONY: ssh-host-keys
ssh-host-keys:
	sudo ssh-keygen -b 1024 -t rsa -f /etc/ssh/ssh_host_key
	sudo ssh-keygen -b 1024 -t rsa -f /etc/ssh/ssh_host_rsa_key 
	sudo ssh-keygen -b 1024 -t dsa -f /etc/ssh/ssh_host_dsa_key

.PHONY: sshd-start
sshd-start:
	sudo service ssh start

.PHONY: ngrok-ssh-start
ngrok-ssh-start:
	ngrok tcp 22

.PHONY: add-sudoers
add-sudoers:
	echo bob     ALL=(ALL:ALL) NOPASSWD:ALL | sudo tee -a /etc/sudoers

.PHONY: file-server
file-server:
	deno run --allow-net --allow-read https://deno.land/std/http/file_server.ts

.PHONY: asdf-sbcl
asdf-sbcl:
	# asdf list
	# asdf plugin add sbcl
	asdf install sbcl 2.0.4
	# asdf current
	asdf global sbcl 2.0.4

.PHONY: get-rush
get-rush:
	# gnu parallel alternative
	go get -u github.com/shenwei356/rush/

.PHONY: will-cite-parallel
will-cite-parallel:
	echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null &

.PHONY: one-click-hugo-run
one-click-hugo-run:
	# yarn
	npx yarn start
	# git push origin master

.PHONY: flask-news-update
flask-news-update:
	git push heroku master
	# git push origin master

.PHONY: heroku-database
heroku-database:
	export DATABASE_URL=`heroku config:get DATABASE_URL -a aaa12w3 -j`

.PHONY: create-venv
create-venv:
	sudo apt install -y python3-venv
	python3 -m venv ~/venv-3.x.x

.PHONY: venv-jupyter
venv-jupyter:
	source ~/venv-3.x.x/bin/activate

.PHONY: create-requirements
create-requirements:
	pip freeze > requirements.txt

.PHONY: clean-requirements
clean-requirements:
	pipreqs --force .

#refresh-just-completion-bash:
#	complete -W "$(just --summary)" just

.PHONY: gpg-encrypt
gpg-encrypt:
	gpg -c filename

.PHONY: gpg-decrypt
gpg-decrypt:
	gpg filename

.PHONY: xtar
xtar:
	[ -f ~/.bashrc ] && cp ~/.bashrc ~/.bashrc.old
	[ -f ~/.bash_profile ] && cp ~/.bash_profile ~/.bash_profile.old
	tar -C ~  -zcvf  x.tar.gz  .ssh .vim .vimrc .tmux.conf .shextra .bashrc  .gitconfig  .bash_profile .env

.PHONY: gobox-sym-encrypt
gobox-sym-encrypt:
	gobox sym-encrypt x.tar.gz x.tar.gz.gobox

.PHONY: gobox-sym-decrypt
gobox-sym-decrypt:
	gobox sym-decrypt x.tar.gz.gobox x.tar.gz

.PHONY: gobox-encrypt
gobox-encrypt: 
	gobox encrypt ~/.ssh/gobox_public ~/.ssh/gobox_private x.tar.gz x.tar.gz.gobox

.PHONY: gobox-decrypt
gobox-decrypt:
	gobox decrypt ~/.ssh/gobox_public ~/.ssh/gobox_private  x.tar.gz.gobox x.tar.gz

.PHONY: gobox-genkey
gobox-genkey:
	gobox genkey ~/.ssh/gobox_public ~/.ssh/gobox_private

.PHONY: gobox
gobox:
	go install github.com/danderson/gobox

.PHONY: xuntar
xuntar:
	cp examples/x.tar.gz.gpg
	gpg x.tar.gz.gpg
	cp -r .ssh .ssh.backup
	tar xf x.tar.gz

.PHONY: common-tools
common-tools: common-base common-tools-build-essentials common-tools-1  common-tools-golang common-tools-fzf common-tools-nodejs common-tools-rust
	#brew install htop fd ripgrep bat tree rush exa procs
	wget https://raw.githubusercontent.com/brushtechnology/fabricate/master/fabricate.py

.PHONY: common-base
common-base:
	sudo apt update
	sudo apt install -y vim tmux openssl miller

.PHONY: common-tools-fzf
common-tools-fzf:
	git clone --depth 1 https://github.com/junegunn/fzf.git ~/.fzf
	~/.fzf/install

.PHONY: common-tools-nodejs
common-tools-nodejs:
	sudo apt install -y nodejs npm
	sudo npm install -g yarn

.PHONY: common-tools-golang
common-tools-golang:
	sudo rm -rf /usr/local/go && cd /usr/local && sudo wget https://golang.org/dl/go1.17.1.linux-amd64.tar.gz &&  sudo tar -C /usr/local -xzf go1.17.1.linux-amd64.tar.gz

.PHONY: common-tools-1
common-tools-1:
	sudo apt install -y  htop   tree  curl wget

.PHONY: common-tools-build-essentials
common-tools-build-essentials:
	sudo apt install -y build-essential

.PHONY: common-tools-rust
common-tools-rust: rustup
	cargo install ripgrep
	#cargo install cargo-edit

.PHONY: rustup
rustup:
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

.PHONY: cargo-edit
cargo-edit:
	cargo install cargo-edit

.PHONY: docker-engine
docker-engine:
	sudo apt-get update
	sudo apt-get install -y apt-transport-https ca-certificates curl gnupg lsb-release
	curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
	echo "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu  `lsb_release -cs` stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
	sudo apt-get update
	sudo apt-get install docker-ce docker-ce-cli containerd.io

.PHONY: docker-post-install
docker-post-install:
	-sudo groupadd docker
	sudo usermod -aG docker $$USER
	newgrp docker
	#docker run hello-world

.PHONY: minikube-install
minikube-install:
	curl -LO https://storage.googleapis.com/minikube/releases/latest/minikube-linux-amd64
	sudo install minikube-linux-amd64 /home/$$USER/bin/minikube
	minikube start
	minikube kubectl -- get po -A
	#minikube  dashboard
	#https://minikube.sigs.k8s.io/docs/start/
	minikube stop

.PHONY: gcloud-install
gcloud-install:
	# ln -s /bin/python3 ~/bin
	wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-359.0.0-linux-x86_64.tar.gz
	tar xf google-cloud-sdk-*.tar.gz
	./google-cloud-sdk/install.sh
	# make sure google-cloud-sdk/bin comes before existing gcloud in /snap/bin or remove /snap/bin/gcloud
	sudo mv /snap/bin/gcloud /snap/bin/gcloud-old
	gcloud init
	gcloud components list
	gcloud components install kubectl
