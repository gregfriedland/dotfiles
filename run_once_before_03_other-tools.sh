#!/bin/sh

echo "### Running: run_once_before_03_other-tools.sh"

printf "%.s=" $(seq 1 $(tput cols))
echo
echo "Installing other tools"

echo "Installing rust"
curl https://sh.rustup.rs -sSf | sh -s -- -y

echo "Installing uv"
curl -LsSf https://astral.sh/uv/install.sh | sh

echo "Installing zoxide"
curl -sSfL https://raw.githubusercontent.com/ajeetdsouza/zoxide/main/install.sh | sh

echo "Installing zellij"
cargo install --locked zellij

echo "Installing starship"
curl -sS https://starship.rs/install.sh | sh -s -- --yes
