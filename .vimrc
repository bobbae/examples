"set viminfo='25,\"50,n~/.viminfo
if has("autocmd")
  au BufReadPost * if line("'\"") > 1 && line("'\"") <= line("$") | exe "normal! g'\"" | endif
endif

set ic
set sm
set smartindent
execute pathogen#infect()
syntax on
filetype plugin indent on
set mouse=a
"colorscheme happy_hacking
"colorscheme jellybean

"let g:Powerline_symbols = 'unicode'
"python from powerline.vim import setup as powerline_setup
"python powerline_setup()
"python del powerline_setup

let g:airline_powerline_fonts = 1
let g:airline#extensions#tabline#enabled = 1
let g:airline_powerline_fonts = 1
if !exists('g:airline_symbols')
  let g:airline_symbols = {}
endif
let g:airline_symbols.space = "\ua0"
"let g:airline_theme='bubblegum'
let g:airline_theme='alduin'
set encoding=utf-8
set t_Co=256
let g:airline#extensions#tmuxline#enabled = 1

set rtp+=~/tabnine-vim
