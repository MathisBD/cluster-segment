import "lib/github.com/athas/vector/vector"
import "lib/github.com/athas/vector/vspace"

module mk_mean_shift (V : vector) (scalar : real) = {
  module VS = mk_vspace V scalar
  type vec = VS.vector

  def argmin [n] (xs : [n]scalar.t) : i64 =
    let (i_min, _) = 
      reduce 
        (\(i, x) (j, y) -> if x scalar.< y then (i, x) else (j, y)) 
        (0, scalar.highest) 
        (zip (iota n) xs)
    in i_min

  -- We use a simple (anisotropic) gaussian kernel.
  let kernel sigma x y = 
    let d2 = VS.((x - y) / sigma) |> VS.quadrance
    in scalar.(exp (neg (from_fraction 1 2 * d2)))

  def MS_step [n] (data : [n]vec) (sigma : vec) (x : vec) : vec =
    -- Compute the weighted average of the points, using the weights given by the kernel.
    let weights = map (kernel sigma x) data
    let avg = 
      map2 (\w y -> VS.scale w y) weights data 
      |> reduce (VS.+) VS.zero
      |> VS.scale scalar.((i32 1) / (sum weights))
    in avg

  -- Compute the mode associated to a point, by iterating MS_step
  def MS_point [n] (data : [n]vec) (sigma : vec) (x : vec) : vec = 
    let (_, mean) = loop (shift, mean) = (scalar.highest, x) while scalar.(shift > (f32 0.01)) do
      let new_mean = MS_step data sigma mean
      let new_shift = VS.((new_mean - mean) / sigma) |> VS.norm 
      in (new_shift, new_mean)
    in mean

  def MS_merge [n] (modes : [n]vec) (sigma : vec) : []vec =
    let enum_modes = zip (iota n) modes
    let (_, merged) = 
      filter (\(i, m) -> 
        -- We discard the mode m if there is another mode m' that is close to m 
        -- and that has a strictly smaller index.
        all (\(i', m') -> 
          let dist = VS.((m - m') / sigma) |> VS.norm 
          in scalar.(dist >= f32 0.5) || i <= i') 
          enum_modes) 
        enum_modes
      |> unzip
    in merged

  def MS [n] (data : [n]vec) (sigma : vec) : ([]vec, [n]i64) = 
    let all_modes = map (MS_point data sigma) data
    -- First merge the modes that are too close.
    let modes = MS_merge all_modes sigma
    -- Then find the best mode for each point.
    let labels = map 
      (\x -> 
        let dists = map (\m -> VS.((x - m) / sigma) |> VS.norm) modes 
        in argmin dists) 
      data
    in (modes, labels)
}

module vector_2 = cat_vector vector_1 vector_1
module vector_3 = cat_vector vector_1 vector_2
module vector_5 = cat_vector vector_3 vector_2

module mean_shift_2D = mk_mean_shift vector_2 f32
entry cluster_2D [n] (data : [n][2]f32) (sigma : [2]f32) : ([][2]f32, [n]i64) =
  let cast 't (x : [2]t) = 
    vector_2.from_array (x :> [vector_2.length]t)
  let (modes, labels) = mean_shift_2D.MS (map cast data) (cast sigma) 
  in (modes |> map vector_2.to_array :> [][2]f32, labels)

module mean_shift_5D = mk_mean_shift vector_5 f32
entry cluster_5D [n] (data : [n][5]f32) (sigma : [5]f32) : ([][5]f32, [n]i64) =
  let cast 't (x : [5]t) = 
    vector_5.from_array (x :> [vector_5.length]t)
  let (modes, labels) = mean_shift_5D.MS (map cast data) (cast sigma) 
  in (modes |> map vector_5.to_array :> [][5]f32, labels)
