Before making a PR, please check that

- [ ] The `./common` library can compile without errors
- [ ] All cases in `./examples` can compile without errors
- [ ] The format script in `format.sh` returns no differnces
    - To run it, you will need python3 installed. Python3 can be loaded on the PACE cluster via `ml python/3.9.12-rkxvr6`. You will also need to add this line to your `.bashrc` file (`export PATH="$PATH:$HOME/.local/bin"`) so PACE can run python3.