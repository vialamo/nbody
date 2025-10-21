# C++ Version

1.  **Prerequisites:** You need a C++ compiler and the SFML 3.0.2 development libraries.
    ```bash
    sudo apt update
    sudo apt install build-essential libsfml-dev
    ```
2.  **Compile:** The `pocketfft` library is header-only and included in the repository.
    ```bash
    g++ nbody.cpp -o nbody -lsfml-graphics -lsfml-window -lsfml-system
    ```
3.  **Run:** Execute the compiled program:
    ```bash
    ./nbody
    ```

Links for manual installation:
* SFML_3.0.2: [`https://www.sfml-dev.org/`](https://www.sfml-dev.org/)
* pocketfft_hdronly.h: [`https://github.com/hayguen/pocketfft`](https://github.com/hayguen/pocketfft)
