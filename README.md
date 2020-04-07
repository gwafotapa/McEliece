# McEliece

A straightforward implementation of the McEliece cryptosystem that I coded to get familiar with Rust.

### Prerequisites

You need Rust and Cargo installed. See:
<https://doc.rust-lang.org/cargo/getting-started/installation.html>

### Installing

Once you have Rust and Cargo installed, all you need to do is compile the binary:

```
cargo build --release
```

This compiles a binary `mceliece` you can run with:

```
target/release/mceliece
```

You can also run it with Cargo:

```
cargo run --release
```

Use `--help` to get a description of the available options:

```
cargo run --release -- --help
```

### Running

You need to supply the binary with a command.
It accepts four commands: `keygen`, `encrypt`, `decrypt` and `plaintext`.

#### `keygen`

Generates a random couple (public key, secret key).
Takes two optional filename arguments to output public and secret keys to.
If none are given, default filenames `public_key.mce` and `secret_key.mce` are used.
Supports two switches:
* `-n LENGTH` sets the Goppa code length (default: 1024).
* `-t CORRECTION_CAPACITY` sets the Goppa code correction capacity (default: 50).

#### plaintext

Generates a random plaintext for the supplied public key.
Takes two optional filename arguments, the public key and the generated plaintext.
If none are given, default filenames `public_key.mce` and `plaintext.mce` are used.

#### encrypt

Encrypts the given plaintext with the supplied public key.
Takes three optional filename arguments for the public key, the plaintext and the generated ciphertext.
If none are given, filenames `public_key.mce`, `plaintext.mce` and `ciphertext.mce` are used.

#### decrypt

Decrypts the given ciphertext with the supplied secret key.
Takes three optional filename arguments for the secret key, the ciphertext and the generated decrypted text.
If none are given, filenames `secret_key.mce`, `ciphertext.mce` and `decrypted.mce` are used.

##### A complete example

To generate public key `pk.mce` and secret key `sk.mce`:

`cargo run --release keygen pk.mce sk.mce`

To generate a random plaintext `m` with the public key `pk.mce`:

`cargo run --release plaintext pk.mce m`

To encrypt plaintext `m` with public key `pk.mce` and output ciphertext `c`:

`cargo run --release encrypt pk.mce m c`

To decrypt ciphertext `c` with secret key `sk.mce` and output the result in `d`:

`cargo run --release decrypt sk.mce c d`

## Running the tests

Just run:

```
cargo test
```

To add debugging information:

```
RUST_LOG=info cargo test
```

## Generating the documentation

Just run:

```
cargo doc
```

Then open `target/doc/mceliece/index.html` in any browser

## Author

**Guillaume Wafo-Tapa** - [gwafotapa](https://github.com/gwafotapa)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Bibliography

* Engelbert, Daniela, Raphael Overbeck, and Arthur Schmidt. "A summary of McEliece-type cryptosystems and their security." Journal of Mathematical Cryptology JMC 1.2 (2007): 151-199.
* Gao, Shuhong, and Daniel Panario. "Tests and constructions of irreducible polynomials over finite fields." Foundations of computational mathematics. Springer, Berlin, Heidelberg, 1997. 346-361.
