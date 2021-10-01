#alias si := system-info

#system-info:
#	@echo from ~/justfile
#	@echo "This is an {{arch()}} machine running {{os()}}".
#	@echo PATH is {{env_var("PATH")}}
#	@echo PATH is $PATH
#	@echo current invocation directory is {{invocation_directory()}}

ssh-host-keys:
	sudo ssh-keygen -b 1024 -t rsa -f /etc/ssh/ssh_host_key
	sudo ssh-keygen -b 1024 -t rsa -f /etc/ssh/ssh_host_rsa_key 
	sudo ssh-keygen -b 1024 -t dsa -f /etc/ssh/ssh_host_dsa_key

sshd-start:
	sudo service ssh start

ngrok-ssh-start:
	ngrok tcp 22

add-sudoers:
	echo bob     ALL=(ALL:ALL) NOPASSWD:ALL | sudo tee -a /etc/sudoers

file-server:
	deno run --allow-net --allow-read https://deno.land/std/http/file_server.ts

asdf-sbcl:
	# asdf list
	# asdf plugin add sbcl
	asdf install sbcl 2.0.4
	# asdf current
	asdf global sbcl 2.0.4

get-rush:
	# gnu parallel alternative
	go get -u github.com/shenwei356/rush/

will-cite-parallel:
	echo 'will cite' | parallel --citation 1> /dev/null 2> /dev/null &

one-click-hugo-run:
	# yarn
	npx yarn start
	# git push origin master

flask-news-update:
	git push heroku master
	# git push origin master

heroku-database:
	export DATABASE_URL=`heroku config:get DATABASE_URL -a aaa12w3 -j`

create-venv:
	sudo apt install -y python3-venv
	python3 -m venv ~/venv-3.x.x

venv-jupyter:
	source ~/venv-3.x.x/bin/activate

create-requirements:
	pip freeze > requirements.txt

clean-requirements:
	pipreqs --force .

#refresh-just-completion-bash:
#	complete -W "$(just --summary)" just

gpg-encrypt:
	gpg -c filename

gpg-decrypt:
	gpg filename

xtar:
	[ -f ~/.bashrc ] && cp ~/.bashrc ~/.bashrc.old
	[ -f ~/.bash_profile ] && cp ~/.bash_profile ~/.bash_profile.old
	tar -C ~  -zcvf  x.tar.gz  .ssh .vim .vimrc .tmux.conf .shextra .bashrc  .gitconfig  .bash_profile .env

gobox-encrypt: 
	gobox encrypt ~/.ssh/gobox_public ~/.ssh/gobox_private x.tar.gz x.tar.gz.gobox

gobox-decrypt:
	gobox decrypt ~/.ssh/gobox_public ~/.ssh/gobox_private  x.tar.gz.gobox x.tar.gz

gobox-genkey:
	gobox genkey ~/.ssh/gobox_public ~/.ssh/gobox_private

gobox:
	go get github.com/danderson/gobox

xuntar:
	cp examples/x.tar.gz.gpg
	gpg x.tar.gz.gpg
	cp -r .ssh .ssh.backup
	tar xf x.tar.gz

common-tools: common-base common-tools-build-essentials common-tools-1  common-tools-golang common-tools-fzf common-tools-nodejs common-tools-rust
	#brew install htop fd ripgrep bat tree rush exa procs
	wget https://raw.githubusercontent.com/brushtechnology/fabricate/master/fabricate.py

common-base:
	sudo apt update
	sudo apt install -y vim tmux openssl miller

common-tools-fzf:
	git clone --depth 1 https://github.com/junegunn/fzf.git ~/.fzf
	~/.fzf/install

common-tools-nodejs:
	sudo apt install -y nodejs npm

common-tools-golang:
	sudo rm -rf /usr/local/go && cd /usr/local && sudo wget https://golang.org/dl/go1.17.1.linux-amd64.tar.gz &&  sudo tar -C /usr/local -xzf go1.17.1.linux-amd64.tar.gz

common-tools-1:
	sudo apt install -y  htop   tree  curl wget

common-tools-build-essentials:
	sudo apt install -y build-essential

common-tools-rust: rustup
	cargo install ripgrep
	#cargo install cargo-edit

rustup:
	curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

cargo-edit:
	cargo install cargo-edit
