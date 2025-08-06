#!/usr/bin/env python3
"""
Compatibility module for Python version differences.

Provides compatibility shims for modules that changed between Python versions.
"""

import html

# Python 3.13 compatibility: cgi module was removed
# Even in earlier versions, cgi.escape was deprecated and moved to html.escape
try:
    import cgi
    # Check if cgi.escape exists (deprecated in 3.8, removed in 3.13)
    if hasattr(cgi, 'escape'):
        escape = cgi.escape
    else:
        escape = html.escape
except ImportError:
    # Python 3.13+: cgi module removed entirely
    import html

    # Create a minimal cgi compatibility shim
    class CGICompat:
        """Minimal compatibility shim for removed cgi module."""

        @staticmethod
        def escape(s, quote=None):
            """HTML escape function compatible with old cgi.escape."""
            if quote is None:
                return html.escape(s)
            else:
                return html.escape(s, quote=quote)

    # Make cgi available as a compatibility layer
    cgi = CGICompat()
    escape = html.escape

# Always use html.escape as the standard approach
escape = html.escape

# Make escape function available at module level
__all__ = ['escape']
