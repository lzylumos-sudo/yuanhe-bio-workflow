"""Stable additive error envelope; legacy error string remains available."""
def envelope(message, error_code='CONFIGURATION_ERROR', run_id=None, stage=None, details=None,
             recoverable=False, suggested_action='Inspect the error and correct the input or environment.'):
    return {'error': str(message), 'error_code': error_code, 'message': str(message),
            'run_id': run_id, 'stage': stage, 'details': details or {},
            'recoverable': recoverable, 'suggested_action': suggested_action}
