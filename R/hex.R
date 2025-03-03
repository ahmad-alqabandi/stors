make_hex <- function() {
  imgurl <- system.file("stors_logo.png", package = "stors")
  hexSticker::sticker(imgurl,
                      s_x = 1.0,
                      s_y = 1.025,
                      s_width = 0.55,
                      s_height = 0.55,
                      package = "",
                      p_x = 0.6,
                      p_y = 1.45,
                      p_color = "#FFFFFF",
                      p_family = "serif",
                      p_fontface = "bold",
                      p_size = 5.0,
                      h_fill = "#FFFFFF",
                      h_color = "#000",
                      url = "ahmad-alqabandi.github.io/stors",
                      u_size = 1.33,
                      u_family = "mono",
                      u_color = "#000",
                      u_x = 0.17,
                      u_y = 0.55,
                      u_angle = -30,
                      white_around_sticker = TRUE,
                      filename = file.path("inst", "stors_hex.png"),
                      dpi = 600L)
  usethis::use_logo(file.path("inst", "stors_hex.png"), geometry = "480x556")
}