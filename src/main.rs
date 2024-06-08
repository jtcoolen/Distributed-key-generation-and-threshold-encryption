//  Ideal class group cryptography
//  Copyright (C) 2024  Julien Coolen
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>

use class_group::{ABDeltaTriple, BinaryQF};
use curv::{arithmetic::Converter, BigInt};

fn main() {
    unsafe {
        class_group::pari_init(10000000, 2);
    }
    let a: BigInt = BigInt::from_str_radix("1347310664179468558147371727982960102805371574927252724399119343247182932538452304549609704350360058405827948976558722087559341859252338031258062288910984654814255199874816496621961922792890687089794760104660404141195904459619180668507135317125790028783030121033883873501532619563806411495141846196437", 10).unwrap();

    let b = BigInt::from(2);
    let delta = -BigInt::from(3) * BigInt::from(201);
    let abd = ABDeltaTriple { a, b, delta };
    let pf = BinaryQF::binary_quadratic_form_disc(&abd);
    let pari_qf = pf.qf_to_pari_qf();
    let pf2 = BinaryQF::pari_qf_to_qf(pari_qf);
    assert_eq!(pf, pf2)
}
