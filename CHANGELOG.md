# Changelog

## [0.9.0](https://github.com/sequinsbio/sequintools/compare/v0.8.5...v0.9.0) (2025-09-09)


### Features

* add BamReader trait ([#233](https://github.com/sequinsbio/sequintools/issues/233)) ([dd7e373](https://github.com/sequinsbio/sequintools/commit/dd7e37339d52712945fc5a6622815e3f57b28fe3))
* add library coverage module ([#235](https://github.com/sequinsbio/sequintools/issues/235)) ([c2645ac](https://github.com/sequinsbio/sequintools/commit/c2645acee593924deeb6dedcb7668ff12a5a1429))
* bedcov command uses new coverage module ([#236](https://github.com/sequinsbio/sequintools/issues/236)) ([d5945e7](https://github.com/sequinsbio/sequintools/commit/d5945e77085eade8f52d931a37b56564eb9d3315))


### Bug Fixes

* bedcov does not error when missing coverage ([#237](https://github.com/sequinsbio/sequintools/issues/237)) ([f5343f3](https://github.com/sequinsbio/sequintools/commit/f5343f359e7449493fecad1676eaa397828ecec1))

## [0.8.5](https://github.com/sequinsbio/sequintools/compare/v0.8.4...v0.8.5) (2025-08-26)


### Bug Fixes

* **deps:** bump anyhow from 1.0.98 to 1.0.99 ([#223](https://github.com/sequinsbio/sequintools/issues/223)) ([0ef481b](https://github.com/sequinsbio/sequintools/commit/0ef481b01969ad87c946684a176ebb0a50e72c75))
* **deps:** bump clap from 4.5.41 to 4.5.43 ([#220](https://github.com/sequinsbio/sequintools/issues/220)) ([5a461e8](https://github.com/sequinsbio/sequintools/commit/5a461e8a3b6e812e317043a04c072d20f85775e4))
* **deps:** bump clap from 4.5.43 to 4.5.45 ([#224](https://github.com/sequinsbio/sequintools/issues/224)) ([6824c1b](https://github.com/sequinsbio/sequintools/commit/6824c1b9a60807672cc4fb2bdbe7b59df103227c))
* **deps:** bump debian from 12.11-slim to 13.0-slim ([#222](https://github.com/sequinsbio/sequintools/issues/222)) ([cccf6d7](https://github.com/sequinsbio/sequintools/commit/cccf6d7c30ac7c98b3cfacfae5b999256428040a))
* **deps:** bump tempfile from 3.20.0 to 3.21.0 ([#225](https://github.com/sequinsbio/sequintools/issues/225)) ([9269c2b](https://github.com/sequinsbio/sequintools/commit/9269c2bf647e549788015a62db2a9171bb52a60f))

## [0.8.4](https://github.com/sequinsbio/sequintools/compare/v0.8.3...v0.8.4) (2025-07-22)


### Bug Fixes

* **deps:** bump clap from 4.5.40 to 4.5.41 ([#214](https://github.com/sequinsbio/sequintools/issues/214)) ([71a3c6a](https://github.com/sequinsbio/sequintools/commit/71a3c6a6c7cb40296adbccc346ce3e66db7f3c6d))
* **deps:** bump rust-htslib from 0.49.0 to 0.50.0 ([#213](https://github.com/sequinsbio/sequintools/issues/213)) ([ebd0932](https://github.com/sequinsbio/sequintools/commit/ebd09328ba0eabe47f2d1ac9e1fbc063f580a377))

## [0.8.3](https://github.com/sequinsbio/sequintools/compare/v0.8.2...v0.8.3) (2025-07-01)


### Bug Fixes

* **deps:** bump clap from 4.5.38 to 4.5.40 ([#211](https://github.com/sequinsbio/sequintools/issues/211)) ([b425620](https://github.com/sequinsbio/sequintools/commit/b425620e46a3ab68bf16aea74f3d24a5d98fe49d))

## [0.8.2](https://github.com/sequinsbio/sequintools/compare/v0.8.1...v0.8.2) (2025-05-28)


### Bug Fixes

* **ci:** run builds on correct OS ([#208](https://github.com/sequinsbio/sequintools/issues/208)) ([58632f6](https://github.com/sequinsbio/sequintools/commit/58632f6354a32de51b49280d88cfad8d9d42922f))

## [0.8.1](https://github.com/sequinsbio/sequintools/compare/v0.8.0...v0.8.1) (2025-05-27)


### Bug Fixes

* **ci:** add musl dependencies to build ([#206](https://github.com/sequinsbio/sequintools/issues/206)) ([dc58325](https://github.com/sequinsbio/sequintools/commit/dc58325bccce040ab50bb1e4f65a42dcd09fcb7e))

## [0.8.0](https://github.com/sequinsbio/sequintools/compare/v0.7.0...v0.8.0) (2025-05-27)


### Features

* add support for CRAM output ([#202](https://github.com/sequinsbio/sequintools/issues/202)) ([9021048](https://github.com/sequinsbio/sequintools/commit/9021048aa066516a7f03211aacfa0b1fe8c051b7)), closes [#199](https://github.com/sequinsbio/sequintools/issues/199)


### Bug Fixes

* abort when summary requested without index ([#201](https://github.com/sequinsbio/sequintools/issues/201)) ([5df5059](https://github.com/sequinsbio/sequintools/commit/5df50590810c867d5fc50f771082fd3a9dd97656))
* reduce crate size to complie with crates.io restrictions ([#200](https://github.com/sequinsbio/sequintools/issues/200)) ([f10d6b1](https://github.com/sequinsbio/sequintools/commit/f10d6b16635b14e5256bc7116625e80355f86613))

## [0.7.0](https://github.com/sequinsbio/sequintools/compare/v0.6.1...v0.7.0) (2025-05-22)


### Features

* allow CRAM input files ([#198](https://github.com/sequinsbio/sequintools/issues/198)) ([53b1a98](https://github.com/sequinsbio/sequintools/commit/53b1a98f969eb387d691fac554978a2f0460febf)), closes [#179](https://github.com/sequinsbio/sequintools/issues/179)
* calibrate by region mean depth with sample bed ([#192](https://github.com/sequinsbio/sequintools/issues/192)) ([e94edc3](https://github.com/sequinsbio/sequintools/commit/e94edc3b088494a94a83b2fbe3bea990858e1c46))
* use multi threads when reading bams ([#194](https://github.com/sequinsbio/sequintools/issues/194)) ([03c9755](https://github.com/sequinsbio/sequintools/commit/03c9755c58ac52f17a0bafb80a79891289b0b5f7))

## [0.6.1](https://github.com/sequinsbio/sequintools/compare/v0.6.0...v0.6.1) (2025-05-13)


### Bug Fixes

* **deps:** bump clap from 4.5.37 to 4.5.38 ([#188](https://github.com/sequinsbio/sequintools/issues/188)) ([129c57c](https://github.com/sequinsbio/sequintools/commit/129c57ce983bf1bf23402134bb1c66c731055ef8))
* **deps:** bump tempfile from 3.19.1 to 3.20.0 ([#189](https://github.com/sequinsbio/sequintools/issues/189)) ([0846711](https://github.com/sequinsbio/sequintools/commit/084671184abc7adf066b84fb1c72cd55e0570754))
* legacy cache api ([#184](https://github.com/sequinsbio/sequintools/issues/184)) ([d393494](https://github.com/sequinsbio/sequintools/commit/d393494327c7be81dc90175fa4ebbb90b2e4c0d4))

## [0.6.0](https://github.com/sequinsbio/sequintools/compare/v0.5.4...v0.6.0) (2025-05-08)


### Features

* adjust calibration algorithm for standard coverage ([#185](https://github.com/sequinsbio/sequintools/issues/185)) ([bdded1c](https://github.com/sequinsbio/sequintools/commit/bdded1ca76bb1ba9ef548c7fff595117f98d7477))


### Bug Fixes

* **deps:** bump anyhow from 1.0.97 to 1.0.98 ([#177](https://github.com/sequinsbio/sequintools/issues/177)) ([2a5cfd0](https://github.com/sequinsbio/sequintools/commit/2a5cfd032b885a482fb0c89b604ad7d4df890b03))
* **deps:** bump clap from 4.5.34 to 4.5.36 ([#178](https://github.com/sequinsbio/sequintools/issues/178)) ([002b534](https://github.com/sequinsbio/sequintools/commit/002b534954f0ce8db1136a8facbc33979d1b59a4))
* **deps:** bump clap from 4.5.36 to 4.5.37 ([#182](https://github.com/sequinsbio/sequintools/issues/182)) ([335850c](https://github.com/sequinsbio/sequintools/commit/335850c0bf36267445abed81c96ad763d1a73fda))
* **deps:** bump rand from 0.9.0 to 0.9.1 ([#181](https://github.com/sequinsbio/sequintools/issues/181)) ([5c095c8](https://github.com/sequinsbio/sequintools/commit/5c095c892fa5729d6ee70c812aaa9f8ee2818383))

## [0.5.4](https://github.com/sequinsbio/sequintools/compare/v0.5.3...v0.5.4) (2025-04-04)


### Bug Fixes

* **deps:** bump clap from 4.5.31 to 4.5.32 ([#171](https://github.com/sequinsbio/sequintools/issues/171)) ([27225c1](https://github.com/sequinsbio/sequintools/commit/27225c18229df2fa2743d72ee6c4e11f18182800))
* **deps:** bump clap from 4.5.32 to 4.5.34 ([#175](https://github.com/sequinsbio/sequintools/issues/175)) ([6b9bf97](https://github.com/sequinsbio/sequintools/commit/6b9bf97e11ba1597f5183f3fc6259ab5ad77c7b8))
* **deps:** bump tempfile from 3.18.0 to 3.19.1 ([#173](https://github.com/sequinsbio/sequintools/issues/173)) ([4313b8d](https://github.com/sequinsbio/sequintools/commit/4313b8de62368859d666e91de9f77a0c177a1925))

## [0.5.3](https://github.com/sequinsbio/sequintools/compare/v0.5.2...v0.5.3) (2025-03-11)


### Bug Fixes

* **deps:** bump anyhow from 1.0.96 to 1.0.97 ([#166](https://github.com/sequinsbio/sequintools/issues/166)) ([6f88668](https://github.com/sequinsbio/sequintools/commit/6f88668efe0d227647789dd661e73521c3191193))
* **deps:** bump clap from 4.5.30 to 4.5.31 ([#165](https://github.com/sequinsbio/sequintools/issues/165)) ([b7af107](https://github.com/sequinsbio/sequintools/commit/b7af1078b4949386c9a145bc6a0f3df45e08ac20))
* **deps:** bump tempfile from 3.17.1 to 3.18.0 ([#170](https://github.com/sequinsbio/sequintools/issues/170)) ([d583d43](https://github.com/sequinsbio/sequintools/commit/d583d4330932a129b01470495776c36c2b545d7d))
* keep sample reads when calibrating ([#169](https://github.com/sequinsbio/sequintools/issues/169)) ([2e3023b](https://github.com/sequinsbio/sequintools/commit/2e3023b1c439c5bff5440e85997fbe0f92f37bf7))

## [0.5.2](https://github.com/sequinsbio/sequintools/compare/v0.5.1...v0.5.2) (2025-02-24)


### Bug Fixes

* **deps:** bump anyhow from 1.0.95 to 1.0.96 ([#162](https://github.com/sequinsbio/sequintools/issues/162)) ([89f60f5](https://github.com/sequinsbio/sequintools/commit/89f60f5376e006c9ad7c77464c14de6b0eb2437c))
* **deps:** bump clap from 4.5.28 to 4.5.29 ([#159](https://github.com/sequinsbio/sequintools/issues/159)) ([fbeb232](https://github.com/sequinsbio/sequintools/commit/fbeb23259072b935a70548e17c5eb4ed8d20613f))
* **deps:** bump clap from 4.5.29 to 4.5.30 ([#161](https://github.com/sequinsbio/sequintools/issues/161)) ([ffa5ef7](https://github.com/sequinsbio/sequintools/commit/ffa5ef7676d4a5b4d2ef88c6ee8f653189607f66))
* **deps:** bump tempfile from 3.16.0 to 3.17.0 ([#158](https://github.com/sequinsbio/sequintools/issues/158)) ([aeaba36](https://github.com/sequinsbio/sequintools/commit/aeaba36311ff1ae4615212caeddbe48af7e127e8))
* **deps:** bump tempfile from 3.17.0 to 3.17.1 ([#163](https://github.com/sequinsbio/sequintools/issues/163)) ([cf841bc](https://github.com/sequinsbio/sequintools/commit/cf841bc7bc6e3b4abf85c05c70d5e1b3d8f65db4))

## [0.5.1](https://github.com/sequinsbio/sequintools/compare/v0.5.0...v0.5.1) (2025-02-10)


### Bug Fixes

* **deps:** bump clap from 4.5.27 to 4.5.28 ([#155](https://github.com/sequinsbio/sequintools/issues/155)) ([af87f1c](https://github.com/sequinsbio/sequintools/commit/af87f1ce6622916129b552873aa32e21ddf26d27))
* **deps:** bump tempfile from 3.14.0 to 3.16.0 ([#149](https://github.com/sequinsbio/sequintools/issues/149)) ([8db866f](https://github.com/sequinsbio/sequintools/commit/8db866f1aeafcbfc2b588bcb290aceaee96da012))
* do not apply a depth limit when calculating mean in calibrate ([#154](https://github.com/sequinsbio/sequintools/issues/154)) ([a6cc8a8](https://github.com/sequinsbio/sequintools/commit/a6cc8a85f92b53ac2ede2d0c1c32f56dbc07fbc7))
* uprev random dependencies ([#152](https://github.com/sequinsbio/sequintools/issues/152)) ([fc3a948](https://github.com/sequinsbio/sequintools/commit/fc3a9489f1f59022b44183ec2c0b04f12f8aa41d))

## [0.5.0](https://github.com/sequinsbio/sequintools/compare/v0.4.0...v0.5.0) (2025-01-31)


### Features

* interpret bedcov max depth 0 as no limit ([#147](https://github.com/sequinsbio/sequintools/issues/147)) ([ddbe085](https://github.com/sequinsbio/sequintools/commit/ddbe085f45cb3181d33d93b6735ee0030682ef9d))

## [0.4.0](https://github.com/sequinsbio/sequintools/compare/v0.3.3...v0.4.0) (2025-01-30)


### Features

* allow setting max depth for bedcov ([#140](https://github.com/sequinsbio/sequintools/issues/140)) ([092e602](https://github.com/sequinsbio/sequintools/commit/092e602fc1b7227dd0396395d055ccf5241ef45c)), closes [#139](https://github.com/sequinsbio/sequintools/issues/139)

## [0.3.3](https://github.com/sequinsbio/sequintools/compare/v0.3.2...v0.3.3) (2025-01-21)


### Features

* allow docker image to be used by nextflow ([#130](https://github.com/sequinsbio/sequintools/issues/130)) ([2b2ab12](https://github.com/sequinsbio/sequintools/commit/2b2ab12017f456b176d98abb878843e89c3f620d)), closes [#129](https://github.com/sequinsbio/sequintools/issues/129)


### Bug Fixes

* **deps:** bump clap from 4.5.23 to 4.5.26 ([#127](https://github.com/sequinsbio/sequintools/issues/127)) ([8711063](https://github.com/sequinsbio/sequintools/commit/8711063d708dd2ed4aef3b2d85774e37d7720f40))

## [0.3.2](https://github.com/sequinsbio/sequintools/compare/v0.3.1...v0.3.2) (2024-12-24)


### Bug Fixes

* do not include component in the release tag ([#122](https://github.com/sequinsbio/sequintools/issues/122)) ([d1a2bd0](https://github.com/sequinsbio/sequintools/commit/d1a2bd05eb3b2e0b5412ee667a9690decfca0cb1))

## [0.3.1](https://github.com/sequinsbio/sequintools/compare/sequintools-v0.3.0...sequintools-v0.3.1) (2024-12-24)


### Bug Fixes

* minimize IO operations within loop ([#114](https://github.com/sequinsbio/sequintools/issues/114)) ([810e2b6](https://github.com/sequinsbio/sequintools/commit/810e2b6fae347144b303b7dd22b742a7d08c9f5d))

## Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
