<p align="center"><img alt="LA-vector" src="http://pages.di.unipi.it/vinciguerra/img/la_vector.svg" height="150"></p>

<p align="center">
LA-vector is a compressed bitvector/container supporting efficient random access and rank queries. It uses novel ways of compressing and accessing data by learning and adapting to data regularities, as described in <a href="http://pages.di.unipi.it/vinciguerra/publication/learned-rank-select/">this research paper</a>.
</p>

<p align="center">
  <a href="https://github.com/gvinciguerra/la_vector/actions?query=workflow%3Abuild"><img src="https://img.shields.io/github/workflow/status/gvinciguerra/la_vector/build" alt="GitHub Workflow Status"></a>
  <a href="https://github.com/gvinciguerra/la_vector/blob/master/LICENSE"><img src="https://img.shields.io/github/license/gvinciguerra/la_vector" alt="License"></a>
  <a href="https://github.com/gvinciguerra/la_vector/stargazers"><img src="https://img.shields.io/github/stars/gvinciguerra/la_vector" alt="GitHub stars"></a>
  <a href="https://github.com/gvinciguerra/la_vector/network/members"><img alt="GitHub forks" src="https://img.shields.io/github/forks/gvinciguerra/la_vector"></a>
</p>

This repo provides a reference C++17 implementation of LA-vector with the following features:

- A container interface via the `operator[]` and `lower_bound` methods.
- A succinct bitvector interface via the `select` and `rank` methods.
- Navigation through iterators.
- Serialisation capabilities.

## Usage

LA-vector is implemented in the two header files inside the `include` directory, and it uses the [sdsl](https://github.com/simongog/sdsl-lite) library.

To compile sdsl, tests, and the [example.cpp](example.cpp) file, use the following commands:

```sh
git clone --recurse-submodules https://github.com/gvinciguerra/la_vector.git
cd la_vector
cmake . -DCMAKE_BUILD_TYPE=Release
make -j8
```

## Contribute

Contributions are welcome. Some ideas:

1. Using vector instructions (e.g. in `la_vector::decode` and `::lower_bound`).
2. Compressing segments and corrections.

## License

This project is licensed under the terms of the Apache License 2.0.

If you use this code for your research, please cite:

> Antonio Boffa, Paolo Ferragina, and Giorgio Vinciguerra. A “learned” approach to quicken and compress rank/select dictionaries. In Proceedings of the SIAM Symposium on Algorithm Engineering and Experiments (ALENEX), 2021.

```bibtex
@inproceedings{Boffa:2021,
  Author = {Boffa, Antonio and Ferragina, Paolo and Vinciguerra, Giorgio},
  Booktitle = {Proceedings of the SIAM Symposium on Algorithm Engineering and Experiments (ALENEX)},
  Title = {A ``learned'' approach to quicken and compress rank/select dictionaries},
  Year = {2021}}
```

The code to reproduce the experiments of the paper is available [here](https://github.com/aboffa/Learned-Rank-Select-ALENEX21).
