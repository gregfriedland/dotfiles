#!/bin/sh

echo "Setting up ZSH as the default shell"
sudo chsh $USER -s /usr/bin/zsh

if [ ! -f /home/greg/.ssh/id_ed25519 ]; then
    echo "Creating SSH key"
    ssh-keygen -t ed25519 -C "greg@rezotx.com"
    echo "Add the following to https://github.com/settings/ssh/new:"
    cat /home/greg/.ssh/id_ed25519.pub
fi
