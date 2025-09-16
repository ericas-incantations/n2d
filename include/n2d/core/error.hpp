#pragma once

#include <string>
#include <utility>

namespace n2d {

enum class ErrorCode {
    None,
    IoError,
    InvalidArgument,
    ParseError,
    Unsupported,
    Internal,
};

struct N2DError {
    ErrorCode code = ErrorCode::Internal;
    std::string message;
};

inline N2DError make_error(ErrorCode code, std::string message) {
    return N2DError{code, std::move(message)};
}

inline N2DError make_io_error(std::string message) {
    return make_error(ErrorCode::IoError, std::move(message));
}

} // namespace n2d

