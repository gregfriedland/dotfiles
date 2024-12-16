#!/bin/sh

echo "### Running: run_once_before_02_configure-zsh.sh"

printf "%.s=" $(seq 1 $(tput cols))
echo
echo "Configuring zsh"

# git clone --depth 1 https://github.com/unixorn/fzf-zsh-plugin.git ${ZSH_CUSTOM:-~/.oh-my-zsh/custom}/plugins/fzf-zsh-plugin
mkdir -p /home/greg/.fzf
touch /home/greg/.fzf/fzf.zsh

printf "%.s=" $(seq 1 $(tput cols))
echo
echo "Setting up zsh"
if [ "$SHELL" != "$(which zsh)" ]; then
    sudo chsh -s $(which zsh) $(whoami)
    echo "Changed default shell to zsh. Changes will take effect after logout."
else
    echo "Shell is already set to zsh"
fi
