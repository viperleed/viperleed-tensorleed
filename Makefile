gcc_fcomp := gfortran
gcc_ccomp := gcc

intel_fcomp := ifort
intel_comp := icc

GCC_FFLAGS := -O3 -fno-finite-math-only
GCC_CFLAGS := -O3 -fno-finite-math-only

INTEL_FFLAGS := -O3
INTEL_CFLAGS := -O3

eeasisss_mod_files := bessl.mod dimen.mod error.mod param.mod poten.mod radeq.mod radii.mod scatt.mod

# Set the sha256 command to be used for checksum verification
SHA256_CMD = sha256sum
ifeq ($(shell uname), Darwin)
        SHA256_CMD = openssl sha256 -r
endif


EEASISS_IMPORTED_ROUTINES_SHA256 := 74ed6b893ac9c6f6125f475ddcca135edb8ddec7373463da44dbcc2afe5f8e35
EEASISS_SHA256 := f1740205de5e4afc5199179514744782b9babe616344f79f7a82e7aeb8696527


intel: eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 correct_sha256
	$(intel_fcomp) eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 -o eeasisss $(INTEL_FFLAGS)
	rm $(eeasisss_mod_files)

gcc: eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 correct_sha256
	$(gcc_fcomp) eeasisss_src/imported_routines.f90 eeasisss_src/eeasisss.f90 -o eeasisss $(GCC_FFLAGS)
	rm $(eeasisss_mod_files)

correct_sha256:
	@IMPORTED_ROUTINES_ACTUAL_SHA256=`cat eeasisss_src/imported_routines.f90 | tr -d '\r' | $(SHA256_CMD)`; \
		case "$$IMPORTED_ROUTINES_ACTUAL_SHA256 " in \
			($(EEASISS_IMPORTED_ROUTINES_SHA256)\ *) : ok ;; \
			(*) echo eeasisss_src/imported_routines.f90 checksum mismatch, expected=\"$(EEASISS_IMPORTED_ROUTINES_SHA256)\" actual=\"$$IMPORTED_ROUTINES_ACTUAL_SHA256\"; \
			exit 1 ;; \
		esac
	@EEASISS_ACTUAL_SHA256=`cat eeasisss_src/eeasisss.f90 | tr -d '\r' | $(SHA256_CMD)`; \
		case "$$EEASISS_ACTUAL_SHA256 " in \
			($(EEASISS_SHA256)\ *) : ok ;; \
			(*) echo eeasisss_src/eeasisss.f90 checksum mismatch, expected=\"$(EEASISS_SHA256)\" actual=\"$$EEASISS_ACTUAL_SHA256\"; \
			exit 1 ;; \
		esac

clean:
	rm eeasisss
	rm $(eeasisss_mod_files)