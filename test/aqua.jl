using Aqua
Aqua.test_all(TypeUtils; stale_deps = !isdefined(Base, :get_extension))
